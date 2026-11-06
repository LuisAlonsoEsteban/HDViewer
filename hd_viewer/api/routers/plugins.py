import dataclasses

from fastapi import APIRouter

from ...core.pipeline import segmentation_service

router = APIRouter()


def describe(plugin):
    return dict(name=plugin.name, title=plugin.title, description=plugin.description, url=plugin.url,
                params=[{key: value for key, value in dataclasses.asdict(param).items()
                         if key != "multiple_subjects"} for param in plugin.params])


@router.get("/plugins")
def plugins():
    """Segmentation methods and their parameter forms. The runner (Docker,
    Apptainer, local) is detected in the background on first use."""
    try:
        service = segmentation_service()
    except Exception as exc:
        return dict(available=False, error=str(exc), runner=None, segmentation=[])
    if service.runner_name is None:
        service.detect_runner(lambda runner: None)
    return dict(available=True, error=None, runner=service.runner_name,
                segmentation=[describe(plugin) for plugin in service.plugins.values()])
