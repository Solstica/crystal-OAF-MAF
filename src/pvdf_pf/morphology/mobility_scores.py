"""Declared spatial hypotheses for where OAF devitrification starts.

Rui et al. constrain the *amount* of RAF/MAF change with temperature but do not
identify which OAF locations become mobile first. These scores are therefore
project hypotheses for sensitivity analysis, not source-derived constitutive laws.

The preferred baseline hypothesis is `iaf_proximal_score`: OAF locations closer
to the isotropic/mobile amorphous region are assigned higher mobility. The reverse
`crystal_proximal_score` is retained as a counterfactual bracket.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import distance_transform_edt

from pvdf_pf.morphology.three_phase import CRYSTAL, OAF, MAF


def _periodic_distance_to(mask: np.ndarray) -> np.ndarray:
    """Approximate periodic Euclidean distance to True pixels by 3x3 tiling."""
    target = np.asarray(mask, dtype=bool)
    if target.ndim != 2:
        raise ValueError("mask must be 2-D")
    if not np.any(target):
        raise ValueError("target mask contains no True pixels")
    nz, nx = target.shape
    tiled = np.tile(target, (3, 3))
    # distance_transform_edt measures distance to zeros; invert target.
    dist = distance_transform_edt(~tiled)
    return dist[nz:2 * nz, nx:2 * nx]


def iaf_proximal_score(three_phase_map: np.ndarray) -> np.ndarray:
    """Higher score for OAF pixels closer to project MAF/source-IAF regions.

    Provenance: PROJECT_SPATIAL_HYPOTHESIS.
    """
    phase = np.asarray(three_phase_map)
    if not set(np.unique(phase)).issubset({CRYSTAL, OAF, MAF}):
        raise ValueError("unsupported phase codes")
    dist = _periodic_distance_to(phase == MAF)
    return -dist


def crystal_proximal_score(three_phase_map: np.ndarray) -> np.ndarray:
    """Counterfactual: higher score for OAF pixels closer to crystal regions."""
    phase = np.asarray(three_phase_map)
    if not set(np.unique(phase)).issubset({CRYSTAL, OAF, MAF}):
        raise ValueError("unsupported phase codes")
    dist = _periodic_distance_to(phase == CRYSTAL)
    return -dist
