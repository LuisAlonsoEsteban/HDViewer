
from .....plugins import runners as _shared
from .....plugins.runners import (RunnerError, THREAD_VARIABLES, DEFAULT_THREADS,
    PROGRESS_PREFIX, CONTAINER_INPUT, CONTAINER_OUTPUT)
from .plugins import SDK_DIR, SEGMENTATION_ROOT


class _Defaults:
    def __init__(self, config, *, sdk_dir=SDK_DIR, segmentation_root=SEGMENTATION_ROOT):
        super().__init__(config, sdk_dir=sdk_dir, segmentation_root=segmentation_root)


class Runner(_Defaults, _shared.Runner):
    pass


class DockerRunner(_shared.DockerRunner, Runner):
    pass


class ApptainerRunner(_shared.ApptainerRunner, Runner):
    pass


class LocalRunner(_shared.LocalRunner, Runner):
    pass


RUNNERS = {cls.name: cls for cls in (DockerRunner, ApptainerRunner, LocalRunner)}


def make_runner(config, *, sdk_dir=SDK_DIR, segmentation_root=SEGMENTATION_ROOT):
    return _shared.make_runner(config, sdk_dir=sdk_dir, segmentation_root=segmentation_root,
                               runner_types=RUNNERS)
