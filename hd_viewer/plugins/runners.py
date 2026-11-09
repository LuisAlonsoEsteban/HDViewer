"""Docker, Apptainer and local runners shared by analysis plugins."""

import json
import os
import shutil
import signal
import subprocess
import sys
import threading
import time
import uuid
from collections import deque
from pathlib import Path

from .jobs import Cancelled, current_job


PROGRESS_PREFIX = "##progress"
CONTAINER_INPUT = "/data/input"
CONTAINER_OUTPUT = "/data/output"

# Numerical libraries start one thread per core by default; with several methods
# running at once on a shared machine that exhausts per-user thread limits
THREAD_VARIABLES = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "NUMEXPR_NUM_THREADS", "TF_NUM_INTRAOP_THREADS", "VHD_SEG_THREADS", "VHD_DECONV_THREADS")
DEFAULT_THREADS = 8
# Image inspection is local; pulls/builds of multi-GB images may take long
INSPECT_TIMEOUT = 60
DEFAULT_IMAGE_TIMEOUT_MINUTES = 60


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
        # Wall-clock limit of one run; a hung method would otherwise hold its slot forever
        self.timeout = _seconds(config.get("timeout_minutes"))

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

    def command(self, plugin, input_dir, output_dir, env, name=None):
        """(argv, environment for the launcher) running the plugin with `env` set inside.
        `name` identifies the run (e.g. the container) so that it can be stopped."""
        raise NotImplementedError

    def stop(self, name):
        """Stop what the process group kill does not reach (e.g. a container)."""

    def redact(self, value):
        for key in ("VHD_CIBERSORTX_USERNAME", "VHD_CIBERSORTX_TOKEN"):
            secret = os.environ.get(key)
            if secret:
                value = value.replace(secret, "[REDACTED]")
        return value

    def run(self, plugin, input_dir, output_dir, on_progress, log_path):
        name = f"vhd-{plugin.name}-{uuid.uuid4().hex[:12]}"
        with open(log_path, "w") as log:
            self.prepare(plugin, on_progress, log)
            argv, env = self.command(plugin, Path(input_dir), Path(output_dir), self.plugin_env(), name=name)
            log.write(self.redact("$ " + " ".join(map(str, argv))) + "\n")
            log.flush()
            on_progress(0, f'Launching {plugin.title} ({self.name}); image download and environment startup may take time')
            tail = self._stream(argv, env, log, on_progress, timeout=self.timeout, stop=lambda: self.stop(name))
        return tail

    def _stream(self, argv, env, log, on_progress=None, *, timeout=None, stop=None):
        """Run argv, logging its output. It is killed (with `stop()`) after `timeout`
        seconds or when the calling thread's job is cancelled."""
        job = current_job()
        if job is not None:
            job.check()
        tail = deque(maxlen=30)
        try:
            # Own process group: a kill reaches apptainer/R/python children too
            proc = subprocess.Popen(
                [str(a) for a in argv], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1, env={**os.environ, **(env or {})}, start_new_session=True,
            )
        except FileNotFoundError as exc:
            raise RunnerError(f"{argv[0]} not found") from exc

        stopped = []
        stopping = threading.Lock()

        def kill(reason):
            with stopping:
                if stopped:
                    return
                stopped.append(reason)
            if stop is not None:
                try:
                    stop()
                except Exception as exc:
                    log.write(f"[runner] could not stop the run: {exc}\n")
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except OSError:
                pass

        timer = threading.Timer(timeout, kill, ("timeout",)) if timeout else None
        if timer is not None:
            timer.daemon = True
            timer.start()
        if job is not None:
            job.on_cancel(lambda: kill("cancelled"))
        try:
            self._read(proc, log, tail, on_progress)
        except BaseException:
            kill("error")
            raise
        finally:
            proc.stdout.close()
            proc.wait()
            if timer is not None:
                timer.cancel()
            if job is not None:
                job.on_cancel(None)
        if stopped:
            if stopped[0] == "cancelled":
                log.write(f"[runner] cancelled: {job.reason}\n")
                log.flush()
                raise Cancelled(f"cancelled: {job.reason}")
            reason = (f"timed out after {timeout / 60:g} minutes and was stopped ({self.name} runner); "
                      "raise timeout_minutes in config.yaml for larger inputs")
            log.write(f"[runner] {reason}\n")
            log.flush()
            raise RunnerError(reason)
        return self._check(proc.returncode, log, tail)

    def _read(self, proc, log, tail, on_progress):
        last_output_update = float('-inf')
        for line in proc.stdout:
            line = self.redact(line)
            log.write(line)
            log.flush()
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
                # Runtime pulls/imports and R startup often emit ordinary output,
                # before the plugin emits its first structured progress message.
                now = time.monotonic()
                if on_progress is not None and now - last_output_update >= 1:
                    on_progress(0, line[:500])
                    last_output_update = now

    def _check(self, code, log, tail):
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

    @property
    def image_timeout(self):
        return _seconds(self.config.get("docker", {}).get("image_timeout_minutes", DEFAULT_IMAGE_TIMEOUT_MINUTES))

    def _inspect(self, image, *options):
        try:
            return subprocess.run([self.binary, "image", "inspect", *options, image],
                                  capture_output=True, text=True, timeout=INSPECT_TIMEOUT)
        except subprocess.TimeoutExpired:
            raise RunnerError(f"docker image inspect {image} timed out: is the Docker daemon responsive?") from None

    def prepare(self, plugin, on_progress, log):
        if self._inspect(plugin.image).returncode == 0:
            return

        # Plugins ship as prebuilt images: pull from the registry named in the manifest
        on_progress(0.0, "Pulling container image (first run only)")
        argv = [self.binary, "pull", plugin.image]
        log.write(self.redact("$ " + " ".join(map(str, argv))) + "\n")
        try:
            self._stream(argv, None, log, timeout=self.image_timeout)
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
            self._stream(argv, None, log, timeout=self.image_timeout)
        except RunnerError as exc:
            raise RunnerError(f"could not pull {plugin.image} ({pull_error}) "
                              f"nor build it from {dockerfile}: {exc}") from None

    def command(self, plugin, input_dir, output_dir, env, name=None):
        argv = [
            self.binary, "run", "--rm", "--cpus", str(self.threads),
            # Results must stay readable/removable by the viewer's user
            "--user", f"{os.getuid()}:{os.getgid()}", "-e", "HOME=/tmp",
            "-v", f"{input_dir.resolve()}:{plugin.input_mount}:ro",
            "-v", f"{output_dir.resolve()}:{plugin.output_mount}",
        ]
        if name:
            argv += ["--name", name]
        memory = self.config.get("memory")
        if memory:
            # Equal swap limit: no swap beyond the RAM cap
            argv += ["--memory", str(memory), "--memory-swap", str(memory)]
        if self.config.get("pids_limit"):
            argv += ["--pids-limit", str(int(self.config["pids_limit"]))]
        for key, value in env.items():
            argv += ["-e", f"{key}={value}"]
        if self.gpu:
            argv += ["--gpus", "all"]
        if not plugin.secret_args:
            return argv + [plugin.image] + list(plugin.command), None
        # Secrets travel in the launcher's environment (`-e NAME` copies it), never in argv;
        # a shell in the image appends them to the image's own entrypoint
        entrypoint = self._entrypoint(plugin.image)
        for variable in plugin.secret_args.values():
            argv += ["-e", variable]
        return (argv + ["--entrypoint", "/bin/sh", plugin.image, "-c", _secret_script(plugin.secret_args),
                        "sh", *entrypoint, *plugin.command],
                _secret_env(plugin.secret_args))

    def _entrypoint(self, image):
        result = self._inspect(image, "--format", "{{json .Config.Entrypoint}}")
        if result.returncode != 0:
            raise RunnerError(f"cannot read the entrypoint of {image}: {result.stderr.strip()}")
        return json.loads(result.stdout or "null") or []

    def stop(self, name):
        # Killing `docker run` leaves the container running: remove it explicitly
        try:
            subprocess.run([self.binary, "rm", "-f", name], capture_output=True, timeout=INSPECT_TIMEOUT)
        except (OSError, subprocess.TimeoutExpired):
            pass


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

    def command(self, plugin, input_dir, output_dir, env, name=None):
        argv = [self.binary, "exec" if plugin.secret_args else "run", "--cleanenv",
                "--bind", f"{input_dir.resolve()}:{plugin.input_mount}:ro,"
                          f"{output_dir.resolve()}:{plugin.output_mount}"]
        if plugin.name in ("spacexr", "music"):
            # --cleanenv still mounts HOME. R otherwise discovers host packages
            # in ~/R, whose compiled dependencies may not exist in the image.
            env = {**env, "R_LIBS_USER": "/nonexistent"}
        for name, value in env.items():
            argv += ["--env", f"{name}={value}"]
        if self.gpu:
            argv.append("--nv")
        # Needs cgroups (v2 with delegation for unprivileged users); off unless configured
        if self.settings.get("memory"):
            argv += ["--memory", str(self.settings["memory"])]
        if not plugin.secret_args:
            return argv + [self.image_uri(plugin)] + list(plugin.command), None
        # APPTAINERENV_<NAME> reaches the container despite --cleanenv, outside argv
        launcher = {f"{prefix}ENV_{name}": value for name, value in _secret_env(plugin.secret_args).items()
                    for prefix in ("APPTAINER", "SINGULARITY")}
        return (argv + [self.image_uri(plugin), "/bin/sh", "-c",
                        _secret_script(plugin.secret_args, "/.singularity.d/runscript"), "sh", *plugin.command],
                launcher)


class LocalRunner(Runner):
    name = "local"

    def available(self):
        return True

    def command(self, plugin, input_dir, output_dir, env, name=None):
        python = (self.config.get("local", {}).get("python") or {}).get(plugin.name) or sys.executable
        pythonpath = os.pathsep.join(filter(None, [str(self.sdk_dir), os.environ.get("PYTHONPATH")]))
        argv = [python, plugin.path / "run.py", "--input", input_dir, "--output", output_dir]
        return argv, {**env, "PYTHONPATH": pythonpath}


def _seconds(minutes):
    return float(minutes) * 60 if minutes else None


def _secret_env(secret_args):
    missing = [name for name in secret_args.values() if not os.environ.get(name)]
    if missing:
        raise RunnerError(f"set {', '.join(missing)} on the server")
    return {name: os.environ[name] for name in secret_args.values()}


def _secret_script(secret_args, program=None):
    """sh -c script running "$@" (after `program`) with flag "$VARIABLE" pairs appended."""
    for name in secret_args.values():
        if not name.isidentifier():
            raise ValueError(f"invalid variable name {name!r}")
    flags = " ".join(f'{flag} "${name}"' for flag, name in secret_args.items())
    return f'exec {program + " " if program else ""}"$@" {flags}'


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
