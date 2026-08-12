"""Audit Rui et al. 2022 main-text Figure 5B against SI Eqs. S2-S4.

The main article calculates x_RAF/x_MAF with crystallinity x_c=0.59, hence their
sum is approximately 0.41. SI Section S2 subsequently rounds eta_cr to 0.60 and
assumes eta_IAF=0.20, eta_OAF~=0.20. Applying SI Eq. S4 literally to the raw
Figure 5B values can therefore violate closure and can give negative eta_MOAF near
Tg (e.g. x_MAF=0.166 at -30 C < eta_IAF=0.20).

This module exposes that inconsistency instead of silently repairing it.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from pvdf_pf.calibration.source_data import load_digitized_curve_csv


@dataclass(frozen=True)
class Figure5BAuditPoint:
    temperature_C: float
    x_raf: float
    x_maf: float
    raw_amorphous_sum: float
    raw_closure_residual_vs_0p41: float
    si_implied_eta_roaf: float
    si_implied_eta_moaf: float
    si_total_fraction: float
    si_closure_residual: float
    si_nonnegative: bool
    literal_si_mapping_valid: bool

    def to_dict(self) -> dict:
        return asdict(self)


def audit_figure5b(
    csv_path: str | Path,
    *,
    sample_state: str = "melt-recrystallized",
    raw_crystallinity: float = 0.59,
    si_crystal_fraction: float = 0.60,
    si_iaf_fraction: float = 0.20,
    raw_closure_tol: float = 0.012,
    si_closure_tol: float = 0.012,
) -> list[Figure5BAuditPoint]:
    points = load_digitized_curve_csv(csv_path)
    sample = sample_state.lower()
    raf = {p.temperature_C: p.value for p in points if p.sample_state.lower() == sample and p.quantity == "x_RAF"}
    maf = {p.temperature_C: p.value for p in points if p.sample_state.lower() == sample and p.quantity == "x_MAF"}
    common = sorted(set(raf) & set(maf))
    if len(common) < 2:
        raise ValueError("Figure 5B requires paired x_RAF/x_MAF source points")

    raw_amorphous = 1.0 - float(raw_crystallinity)
    out: list[Figure5BAuditPoint] = []
    for T in common:
        xr = float(raf[T])
        xm = float(maf[T])
        raw_sum = xr + xm
        eta_roaf = xr
        eta_moaf = xm - float(si_iaf_fraction)
        total = float(si_crystal_fraction) + eta_roaf + eta_moaf + float(si_iaf_fraction)
        nonnegative = eta_roaf >= -1e-12 and eta_moaf >= -1e-12
        raw_ok = abs(raw_sum - raw_amorphous) <= raw_closure_tol
        si_ok = abs(total - 1.0) <= si_closure_tol
        out.append(
            Figure5BAuditPoint(
                temperature_C=float(T),
                x_raf=xr,
                x_maf=xm,
                raw_amorphous_sum=raw_sum,
                raw_closure_residual_vs_0p41=raw_sum - raw_amorphous,
                si_implied_eta_roaf=eta_roaf,
                si_implied_eta_moaf=eta_moaf,
                si_total_fraction=total,
                si_closure_residual=total - 1.0,
                si_nonnegative=bool(nonnegative),
                literal_si_mapping_valid=bool(raw_ok and si_ok and nonnegative),
            )
        )
    return out


def source_anchor_checks(csv_path: str | Path) -> dict[str, bool]:
    """Check values explicitly reported in the main-text paragraph."""
    points = load_digitized_curve_csv(csv_path)
    lookup = {(p.temperature_C, p.quantity): p.value for p in points}
    expected = {
        (-55.0, "x_RAF"): 0.400,
        (-55.0, "x_MAF"): 0.008,
        (-45.2, "x_RAF"): 0.331,
        (-45.2, "x_MAF"): 0.079,
        (-30.0, "x_RAF"): 0.244,
        (-30.0, "x_MAF"): 0.166,
        (40.0, "x_RAF"): 0.014,
        (40.0, "x_MAF"): 0.396,
    }
    return {
        f"{T:g}_{q}": bool((T, q) in lookup and np.isclose(lookup[(T, q)], val, atol=5e-4))
        for (T, q), val in expected.items()
    }
