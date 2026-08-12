"""Spatial ROAF/MOAF allocation without inventing a devitrification mechanism.

The source model constrains *how much* OAF is rigid/mobile at a temperature once
x_RAF(T) and x_MAF(T) are known. It does not state which OAF locations devitrify
first. This module therefore requires an externally supplied `mobility_score` field.
The score may later come from MD, a calibrated local-environment model, or a stated
hypothesis used only for sensitivity analysis.

No default spatial score is provided intentionally.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pvdf_pf.calibration.phase_fraction_state import SourcePhaseFractionState
from pvdf_pf.morphology.three_phase import CRYSTAL, OAF, MAF

SOURCE_CRYSTAL = 0
SOURCE_ROAF = 1
SOURCE_MOAF = 2
SOURCE_IAF = 3
SOURCE_PHASE_NAMES = {
    SOURCE_CRYSTAL: "crystal",
    SOURCE_ROAF: "ROAF",
    SOURCE_MOAF: "MOAF",
    SOURCE_IAF: "IAF",
}


@dataclass(frozen=True)
class SpatialAllocationReport:
    target_mobile_fraction_within_oaf: float
    realized_mobile_fraction_within_oaf: float
    n_oaf_pixels: int
    n_moaf_pixels: int
    score_min_selected: float | None
    score_max_unselected: float | None


def allocate_source_oaf_subphases(
    three_phase_map: np.ndarray,
    *,
    fraction_state: SourcePhaseFractionState,
    mobility_score: np.ndarray,
) -> tuple[np.ndarray, SpatialAllocationReport]:
    """Convert a crystal/OAF/MAF map into crystal/ROAF/MOAF/IAF.

    The existing project MAF pixels are mapped to source IAF *only inside this
    source-model reconstruction*. This does not establish a general physical
    identity between project MAF and Rui-2022 IAF.

    OAF pixels with the largest externally supplied mobility scores are assigned
    to MOAF until the target `mobile_fraction_within_oaf` is reached. Ties are
    broken deterministically by flattened array index.
    """
    phase = np.asarray(three_phase_map)
    score = np.asarray(mobility_score, dtype=float)
    if phase.shape != score.shape:
        raise ValueError("three_phase_map and mobility_score must have identical shapes")
    if np.any(~np.isfinite(score)):
        raise ValueError("mobility_score must be finite everywhere")
    allowed = {CRYSTAL, OAF, MAF}
    if not set(np.unique(phase)).issubset(allowed):
        raise ValueError("three_phase_map contains unsupported phase codes")

    oaf_flat = np.flatnonzero(phase.ravel() == OAF)
    n_oaf = int(oaf_flat.size)
    if n_oaf == 0:
        raise ValueError("three_phase_map contains no OAF pixels")

    target = float(fraction_state.mobile_fraction_within_oaf)
    if not (0.0 <= target <= 1.0):
        raise ValueError("fraction_state has invalid mobile OAF fraction")
    n_mobile = int(np.rint(target * n_oaf))
    n_mobile = min(max(n_mobile, 0), n_oaf)

    out = np.empty(phase.shape, dtype=np.int8)
    out[phase == CRYSTAL] = SOURCE_CRYSTAL
    out[phase == MAF] = SOURCE_IAF
    out[phase == OAF] = SOURCE_ROAF

    selected_score_min = None
    unselected_score_max = None
    if n_mobile > 0:
        score_flat = score.ravel()
        # Lexsort: primary key descending score, secondary key ascending flat index.
        order = np.lexsort((oaf_flat, -score_flat[oaf_flat]))
        selected = oaf_flat[order[:n_mobile]]
        out.ravel()[selected] = SOURCE_MOAF
        selected_score_min = float(np.min(score_flat[selected]))
        if n_mobile < n_oaf:
            remaining = oaf_flat[order[n_mobile:]]
            unselected_score_max = float(np.max(score_flat[remaining]))
    else:
        unselected_score_max = float(np.max(score.ravel()[oaf_flat]))

    realized = float(n_mobile / n_oaf)
    report = SpatialAllocationReport(
        target_mobile_fraction_within_oaf=target,
        realized_mobile_fraction_within_oaf=realized,
        n_oaf_pixels=n_oaf,
        n_moaf_pixels=n_mobile,
        score_min_selected=selected_score_min,
        score_max_unselected=unselected_score_max,
    )
    return out, report


def source_phase_fractions(source_map: np.ndarray) -> dict[str, float]:
    """Return realized fractions in a source-model four-component map."""
    phase = np.asarray(source_map)
    n = phase.size
    if n == 0:
        raise ValueError("source_map must not be empty")
    return {
        SOURCE_PHASE_NAMES[code].lower(): float(np.count_nonzero(phase == code) / n)
        for code in SOURCE_PHASE_NAMES
    }
