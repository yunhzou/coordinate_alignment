"""Reusable augmented fragment detection."""

from .detection import detect_fragments
from .models import (
    FragmentCandidate,
    FragmentDetectionConfig,
    FragmentDetectionResult,
)

__all__ = [
    "FragmentCandidate",
    "FragmentDetectionConfig",
    "FragmentDetectionResult",
    "detect_fragments",
]
