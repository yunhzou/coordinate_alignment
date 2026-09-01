"""Assembly of detected fragment candidates into target covers."""

from .coverage import assemble_fragment_cover
from .models import RetroAssembly, RetroAssemblySearchResult
from .ownership import OwnershipResolution, resolve_overlapping_ownership

__all__ = [
    "RetroAssembly",
    "RetroAssemblySearchResult",
    "OwnershipResolution",
    "assemble_fragment_cover",
    "resolve_overlapping_ownership",
]
