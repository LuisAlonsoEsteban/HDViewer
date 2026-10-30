"""Docker, Apptainer and local runners shared by analysis plugins."""

import os
import shutil
import subprocess
import sys
from collections import deque
from pathlib import Path


PROGRESS_PREFIX = "##progress"
CONTAINER_INPUT = "/data/input"
CONTAINER_OUTPUT = "/data/output"

# Numerical libraries start one thread per core by default; with several methods
# running at once on a shared machine that exhausts per-user thread limits
THREAD_VARIABLES = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "NUMEXPR_NUM_THREADS", "TF_NUM_INTRAOP_THREADS", "VHD_SEG_THREADS", "VHD_DECONV_THREADS")
DEFAULT_THREADS = 8


class RunnerError(RuntimeError):
    pass


class Runner:
    name = None

    def __init__(self, config, *, sdk_dir, segmentation_root):
        self.sdk_dir = Path(sdk_dir)
        self.segmentation_root = Path(segmentation_root)
        self.config = config
        self.gpu = bool(config.get("gpu", False))
        threads = config.get("threads_per_run", DEFAULT_THREADS)
        self.threads = max(1, int(threads if threads != "auto" else os.cpu_count() or 1))

    def plugin_env(self):
        """Environment set inside every plugin run."""
        env = {name: str(self.threads) for name in THREAD_VARIABLES}
        env["TF_NUM_INTEROP_THREADS"] = "2"
        env["VHD_DECONV_GPU"] = "1" if self.gpu else "0"
        # Ignore pip --user packages (~/.local): Apptainer mounts the home directory,
        # and a host numpy 2 there shadowed the image's numpy 1 (TensorFlow 2.15)
        env["PYTHONNOUSERSITE"] = "1"
        return env

    def available(self):
        raise NotImplementedError

    def prepare(self, plugin, on_progress, log):
        """Make the plugin runnable (e.g. pull its image). Called before each run."""

    def command(self, plugin, input_dir, output_dir, env):
        """(argv, environment for the launcher) running the plugin with `env` set inside."""
        raise NotImplementedError

    def redact(self, value):
        for key in ("VHD_CIBERSORTX_USERNAME", "VHD_CIBERSORTX_TOKEN"):
            secret = os.environ.get(key)
            if secret:
                value = value.replace(secret, "[REDACTED]")
        return value

    def run(self, plugin, input_dir, output_dir, on_progress, log_path):
        with open(log_path, "w") as log:
            self.prepare(plugin, on_progress, log)
            argv, env = self.command(plugin, Path(input_dir), Path(output_dir), self.plugin_env())
            log.write(self.redact("$ " + " ".join(map(str, argv))) + "\n")
            log.flush()
            tail = self._stream(argv, env, log, on_progress)
        return tail

    def _stream(self, argv, env, log, on_progress=None):
        tail = deque(maxlen=30)
        try:
            proc = subprocess.Popen(
                [str(a) for a in argv], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1, env={**os.environ, **(env or {})},
            )
        except FileNotFoundError as exc:
            raise RunnerError(f"{argv[0]} not found") from exc

        with proc:
            for line in proc.stdout:
                line = self.redact(line)
                log.write(line)
                line = line.rstrip()
                if line.startswith(PROGRESS_PREFIX):
                    parts = line.split(maxsplit=2)
                    try:
                        fraction = float(parts[1])
                    except (IndexError, ValueError):
                        continue
                    if on_progress is not None:
                        on_progress(fraction, parts[2] if len(parts) > 2 else "")
                elif line:
                    tail.append(line)
        code = proc.returncode
        if code != 0:
            reason = f"exited with code {code} ({self.name} runner)"
            # SIGKILL: -9 for a direct child, 128 + 9 through apptainer/docker.
            # The kernel kills a process over its memory limit this way, so it
            # cannot print anything itself.
            if code in (-9, 137):
                reason += ": killed, most likely out of memory (check the RAM limit of the session or job)"
            # The log otherwise just stops: say why it ended
            log.write(f"[runner] {reason}\n")
            log.flush()
            raise RunnerError(reason + ":\n" + "\n".join(tail))
        log.flush()
        return list(tail)


class DockerRunner(Runner):
    name = "docker"

    @property
    def binary(self):
        return self.config.get("docker", {}).get("binary", 'docker')

    def available(self):
        if shutil.which(self.binary) is None:
            return False
        try:
            return subprocess.run([self.binary, "info"], capture_output=True, timeout=15).returncode == 0
        except (OSError, subprocess.TimeoutExpired):
            return False

    def prepare(self, plugin, on_progress, log):
        exists = subprocess.run([self.binary, "image", "inspect", plugin.image],
                                capture_output=True).returncode == 0
        if exists:
            return

        # Plugins ship as prebuilt images: pull from the registry named in the manifest
        on_progress(0.0, "Pulling container image (first run only)")
        argv = [self.binary, "pull", plugin.image]
        log.write(self.redact("$ " + " ".join(map(str, argv))) + "\n")
        try:
            self._stream(argv, None, log)
            return
        except RunnerError as exc:
            pull_error = exc

        # Unpublished plugins (e.g. under development) can still be built from source
        dockerfile = plugin.path / "Dockerfile"
        if not dockerfile.exists():
            raise RunnerError(f"could not pull {plugin.image}: {pull_error}")
        on_progress(0.0, "Image not in a registry: building it (first run only)")
        argv = [self.binary, "build", "-t", plugin.image, "-f", dockerfile, self.segmentation_root]
        log.write(self.redact("$ " + " ".join(map(str, argv))) + "\n")
        try:
            self._stream(argv, None, log)
        except RunnerError as exc:
            raise RunnerError(f"could not pull {plugin.image} ({pull_error}) "
                              f"nor build it from {dockerfile}: {exc}") from None

    def command(self, plugin, input_dir, output_dir, env):
        argv = [
            self.binary, "run", "--rm", "--cpus", str(self.threads),
            # Results must stay readable/removable by the viewer's user
            "--user", f"{os.getuid()}:{os.getgid()}", "-e", "HOME=/tmp",
            "-v", f"{input_dir.resolve()}:{plugin.input_mount}:ro",
            "-v", f"{output_dir.resolve()}:{plugin.output_mount}",
        ]
        for name, value in env.items():
            argv += ["-e", f"{name}={value}"]
        if self.gpu:
            argv += ["--gpus", "all"]
        return argv + [plugin.image] + list(plugin.command), None


class ApptainerRunner(Runner):
    name = "apptainer"

    @property
    def settings(self):
        return self.config.get("apptainer", {})

    @property
    def binary(self):
        for candidate in (self.settings.get("binary", 'apptainer'), "apptainer", "singularity"):
            if shutil.which(candidate):
                return candidate
        return None

    def available(self):
        return self.binary is not None

    def image_uri(self, plugin):
        configured = (self.settings.get("images") or {}).get(plugin.name)
        if configured:
            return configured
        images_dir = self.segmentation_root / self.settings.get("images_dir", 'images')
        sif = images_dir / f"{plugin.name}.sif"
        if sif.exists():
            return str(sif)
        # Apptainer pulls the registry image and caches it as a SIF on first use
        return f"docker://{plugin.image}"

    def command(self, plugin, input_dir, output_dir, env):
        argv = [self.binary, "run", "--cleanenv",
                "--bind", f"{input_dir.resolve()}:{plugin.input_mount}:ro,"
                          f"{output_dir.resolve()}:{plugin.output_mount}"]
        for name, value in env.items():
            argv += ["--env", f"{name}={value}"]
        if self.gpu:
            argv.append("--nv")
        return argv + [self.image_uri(plugin)] + list(plugin.command), None


class LocalRunner(Runner):
    name = "local"

    def available(self):
        return True

    def command(self, plugin, input_dir, output_dir, env):
        python = (self.config.get("local", {}).get("python") or {}).get(plugin.name) or sys.executable
        pythonpath = os.pathsep.join(filter(None, [str(self.sdk_dir), os.environ.get("PYTHONPATH")]))
        argv = [python, plugin.path / "run.py", "--input", input_dir, "--output", output_dir]
        return argv, {**env, "PYTHONPATH": pythonpath}


RUNNERS = {cls.name: cls for cls in (DockerRunner, ApptainerRunner, LocalRunner)}


def make_runner(config, *, sdk_dir, segmentation_root, env_var="VHD_SEG_RUNNER", kind="segmentation", runner_types=None):
    runners = runner_types or RUNNERS
    choice = os.environ.get(env_var) or config.get("runner", 'auto')
    if choice != "auto":
        if choice not in runners:
            raise ValueError(f"unknown {kind} runner '{choice}'; use auto or {sorted(runners)}")
        return runners[choice](config, sdk_dir=sdk_dir, segmentation_root=segmentation_root)
    for cls in (runners["docker"], runners["apptainer"]):
        runner = cls(config, sdk_dir=sdk_dir, segmentation_root=segmentation_root)
        if runner.available():
            return runner
    return runners["local"](config, sdk_dir=sdk_dir, segmentation_root=segmentation_root)
