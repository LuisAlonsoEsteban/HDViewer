"""Pluggable spatial deconvolution without viewer dependencies."""
from .plugins import Param, Plugin, discover
from .service import DeconvolutionService, DeconvolutionResult, validate_result
from .api.sdk.nnls import deconvolve, _normalise

__all__ = ['Param', 'Plugin', 'discover', 'DeconvolutionService', 'DeconvolutionResult', 'validate_result', 'deconvolve']
