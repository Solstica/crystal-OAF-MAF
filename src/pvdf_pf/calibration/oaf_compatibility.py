"""Compatibility audit for a closed same-state OAF response decomposition.

v0.1.10 provides an effective OAF response required to close the measured
BOPVDF film dielectric constant under the source-anchored structural proxy

    f_crystal = 0.52, f_OAF = 0.28, f_IAF = 0.20.

Rui-2022 Figure 5B provides a temperature-dependent RAF/MAF devitrification
trend, but for a different melt-recrystallized donor state.  v0.1.6 therefore
uses that curve only to define a normalized progress q(T) in [0, 1].

Rui-2022 Figure S2 provides an inverse epsilon_MOAF(T), again through the
cross-state SI construction audited in v0.1.9.

This module asks whether these two transferred quantities can satisfy the
closed same-state mixture identity

    epsilon_OAF,eff = (1-q) epsilon_ROAF + q epsilon_MOAF.

The audit deliberately does not force a fit.  It exposes the epsilon_ROAF that
would be required at q<1, the endpoint incompatibility at q=1, and a residual
under a constant-ROAF anchor hypothesis.  None of these diagnostics is promoted
to an intrinsic material constant or TDGL background permittivity.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from pvdf_pf.calibration.project_devitrification import build_regularized_oaf_state
from pvdf_pf.calibration.source_state_transfer import closed_oaf_effective_response


@dataclass(frozen=True)
class OAFCompatibilityPoint:
    sample_state: str
    temperature_C: float
    q_donor_progress: float
    epsilon_oaf_effective_closed: float
    epsilon_moaf_figure_s2: float
    epsilon_roaf_required: float | None
    endpoint_gap_if_q1: float | None
    epsilon_roaf_anchor_lowT: float
    epsilon_oaf_pred_constant_roaf_anchor: float
    missing_effective_response_constant_roaf_anchor: float
    q_required_constant_roaf_anchor: float | None
    convex_mixture_compatible_with_anchor: bool | None
    provenance: str = "CROSS_STATE_TRANSFER_COMPATIBILITY_AUDIT"

    def to_dict(self) -> dict:
        return asdict(self)


def audit_oaf_compatibility(
    rui2022_dir: str | Path,
    *,
    sample_state: str,
    temperature_C: float,
    low_anchor_C: float = -30.0,
    q_one_tol: float = 1e-10,
) -> OAFCompatibilityPoint:
    root = Path(rui2022_dir)
    sample = str(sample_state).strip().lower()
    T = float(temperature_C)

    target = closed_oaf_effective_response(root, sample_state=sample, temperature_C=T)
    anchor = closed_oaf_effective_response(root, sample_state=sample, temperature_C=low_anchor_C)
    donor = build_regularized_oaf_state(
        root / "figure_5b_maintext_digitized.csv",
        temperature_C=T,
    )

    q = float(donor.q_devitrification)
    eps_eff = float(target.epsilon_oaf_effective_required)
    eps_m = float(target.epsilon_moaf_figure_s2)
    eps_r_anchor = float(anchor.epsilon_oaf_effective_required)

    eps_roaf_required = None
    endpoint_gap = None
    if q < 1.0 - q_one_tol:
        eps_roaf_required = float((eps_eff - q * eps_m) / (1.0 - q))
    else:
        # At q=1 the ROAF coefficient vanishes, so no ROAF value can repair a
        # mismatch between epsilon_OAF,eff and epsilon_MOAF.
        endpoint_gap = float(eps_eff - eps_m)

    eps_pred_anchor = float((1.0 - q) * eps_r_anchor + q * eps_m)
    missing = float(eps_eff - eps_pred_anchor)

    q_required = None
    convex_ok = None
    denom = eps_m - eps_r_anchor
    if abs(denom) > 1e-12:
        q_required = float((eps_eff - eps_r_anchor) / denom)
        convex_ok = bool(-1e-10 <= q_required <= 1.0 + 1e-10)

    return OAFCompatibilityPoint(
        sample_state=sample,
        temperature_C=T,
        q_donor_progress=q,
        epsilon_oaf_effective_closed=eps_eff,
        epsilon_moaf_figure_s2=eps_m,
        epsilon_roaf_required=eps_roaf_required,
        endpoint_gap_if_q1=endpoint_gap,
        epsilon_roaf_anchor_lowT=eps_r_anchor,
        epsilon_oaf_pred_constant_roaf_anchor=eps_pred_anchor,
        missing_effective_response_constant_roaf_anchor=missing,
        q_required_constant_roaf_anchor=q_required,
        convex_mixture_compatible_with_anchor=convex_ok,
    )


def compatibility_series(
    rui2022_dir: str | Path,
    *,
    sample_state: str,
    temperatures_C: tuple[float, ...] = (-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0),
) -> list[OAFCompatibilityPoint]:
    return [
        audit_oaf_compatibility(
            rui2022_dir,
            sample_state=sample_state,
            temperature_C=T,
        )
        for T in temperatures_C
    ]


def summarize_compatibility(points: list[OAFCompatibilityPoint]) -> dict:
    if not points:
        raise ValueError("points must not be empty")

    missing = np.asarray([p.missing_effective_response_constant_roaf_anchor for p in points], dtype=float)
    q_donor = np.asarray([p.q_donor_progress for p in points], dtype=float)
    q_required = np.asarray([
        np.nan if p.q_required_constant_roaf_anchor is None else p.q_required_constant_roaf_anchor
        for p in points
    ], dtype=float)
    roaf_required = np.asarray([
        np.nan if p.epsilon_roaf_required is None else p.epsilon_roaf_required
        for p in points
    ], dtype=float)
    endpoint_gaps = [p.endpoint_gap_if_q1 for p in points if p.endpoint_gap_if_q1 is not None]
    incompatible_anchor = sum(
        1 for p in points if p.convex_mixture_compatible_with_anchor is False
    )

    return {
        "n_points": len(points),
        "q_donor_min": float(np.min(q_donor)),
        "q_donor_max": float(np.max(q_donor)),
        "required_roaf_min_where_identifiable": float(np.nanmin(roaf_required)),
        "required_roaf_max_where_identifiable": float(np.nanmax(roaf_required)),
        "max_abs_missing_response_constant_roaf_anchor": float(np.max(np.abs(missing))),
        "mean_abs_missing_response_constant_roaf_anchor": float(np.mean(np.abs(missing))),
        "n_constant_anchor_convex_incompatible": int(incompatible_anchor),
        "q_required_min_finite": float(np.nanmin(q_required)),
        "q_required_max_finite": float(np.nanmax(q_required)),
        "endpoint_q1_gaps": [float(x) for x in endpoint_gaps],
        "status": (
            "COMPATIBLE" if incompatible_anchor == 0 and all(abs(x) < 1e-6 for x in endpoint_gaps)
            else "TRANSFERRED_Q_AND_FIGURE_S2_MOAF_NOT_JOINTLY_COMPATIBLE"
        ),
    }
