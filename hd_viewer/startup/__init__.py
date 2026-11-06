"""Data loading: readers, samples and operations (Bokeh-free subpackages), plus
the Bokeh loading screen in `controller`.

The controller is imported lazily, so the backend subpackages (sample, data,
operations, utils) import without Bokeh.
"""
import importlib
import pkgutil

_SUBMODULES = {module.name for module in pkgutil.iter_modules(__path__)}


def __getattr__(name):
    if name in _SUBMODULES:
        return importlib.import_module(f"{__name__}.{name}")
    if name.startswith("__"):
        raise AttributeError(name)
    controller = importlib.import_module(f"{__name__}.controller")
    try:
        return getattr(controller, name)
    except AttributeError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
