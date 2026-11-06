import bokeh
from fastapi import APIRouter

from ... import settings
from ...core.branding import BRANDING
from ...startup.sample._sample import Sample

router = APIRouter()


@router.get("/config")
def config():
    return dict(
        branding=BRANDING,
        bokeh_version=bokeh.__version__,
        tutorial_url="/static/assets/tutorial.html",
        use_case_url="/static/assets/biological-use-case.html",
        upload=dict(max_bytes=settings.UPLOAD_MAX_BYTES, chunk_bytes=settings.UPLOAD_CHUNK_BYTES,
                    bin_sizes=list(settings.UPLOAD_BIN_SIZES), default_bin_size=settings.UPLOAD_DEFAULT_BIN_SIZE),
        expression_scales=[dict(value=key, label=label) for key, label in Sample.EXPRESSION_SCALES.items()],
        roi=dict(default_radius=settings.ROI_DEFAULT_RADIUS, preview_debounce_ms=settings.ROI_PREVIEW_DEBOUNCE_MS),
    )
