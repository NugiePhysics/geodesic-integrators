"""Spacetime metrics as compiled kernels (index conventions in :mod:`geoint.metrics.base`)."""

from .base import Metric, MetricKernels, StopSurface
from .schwarzschild import Schwarzschild

__all__ = ["Metric", "MetricKernels", "Schwarzschild", "StopSurface"]
