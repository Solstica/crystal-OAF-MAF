"""Literature-constrained calibration utilities for PVDF phase-field models."""

from .targets import BOPVDF_2021_TARGETS, YANG_2015_TARGETS, MACROMOLECULES_2022_ASSUMPTIONS
from .dielectric import effective_permittivity, infer_unknown_phase_permittivity

__all__ = [
    "BOPVDF_2021_TARGETS",
    "YANG_2015_TARGETS",
    "MACROMOLECULES_2022_ASSUMPTIONS",
    "effective_permittivity",
    "infer_unknown_phase_permittivity",
]
