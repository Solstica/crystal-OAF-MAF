"""Project-regularized OAF devitrification state driven by Rui 2022 Figure 5B.

Why a regularization is needed
-----------------------------
The main-text Figure 5B is calorimetric and uses x_c=0.59, so x_RAF+x_MAF~=0.41.
SI Section S2 later uses rounded structural fractions eta_cr=0.60, eta_IAF=0.20,
eta_OAF~=0.20 and writes x_RAF=eta_ROAF, x_MAF=eta_MOAF+eta_IAF above Tg.
Literal substitution is internally inconsistent near Tg (e.g. at -30 C,
x_MAF=0.166<eta_IAF=0.20) and has a ~0.01 closure mismatch at high T.

The project therefore uses Figure 5B only as a *devitrification progress constraint*.
Between -30 and 40 C, the loss of RAF and gain of MAF are normalized to a progress
q(T) in [0,1]. The SI structural OAF total (0.20) is then partitioned as

    eta_ROAF = eta_OAF * (1-q)
    eta_MOAF = eta_OAF * q.

This is PROJECT_REGULARIZATION_HYPOTHESIS, not an equation reported by Rui et al.
It preserves the measured temperature trend without negative phase fractions.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from pvdf_pf.calibration.phase_fraction_state import SourcePhaseFractionState
from pvdf_pf.calibration.source_data import load_digitized_curve_csv


@dataclass(frozen=True)
class RegularizedOAFDevitrificationState:
    temperature_C: float
    x_raf_raw: float
    x_maf_raw: float
    q_from_raf_loss: float
    q_from_maf_gain: float
    q_devitrification: float
    q_internal_disagreement: float
    crystal: float
    roaf: float
    moaf: float
    iaf: float
    provenance: str = "PROJECT_REGULARIZATION_HYPOTHESIS"

    @property
    def total_oaf(self) -> float:
        return float(self.roaf + self.moaf)

    def to_dict(self) -> dict:
        return asdict(self)

    def as_source_fraction_state(self, *, sample_state: str) -> SourcePhaseFractionState:
        return SourcePhaseFractionState(
            sample_state=str(sample_state).strip().lower(),
            temperature_C=self.temperature_C,
            x_raf=self.x_raf_raw,
            x_maf=self.x_maf_raw,
            crystal=self.crystal,
            roaf=self.roaf,
            moaf=self.moaf,
            iaf=self.iaf,
            mobile_fraction_within_oaf=self.q_devitrification,
        )


def _paired_curves(csv_path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    points = load_digitized_curve_csv(csv_path)
    sample = "melt-recrystallized"
    raf = {float(p.temperature_C): float(p.value) for p in points if p.sample_state.lower() == sample and p.quantity == "x_RAF"}
    maf = {float(p.temperature_C): float(p.value) for p in points if p.sample_state.lower() == sample and p.quantity == "x_MAF"}
    common = sorted(set(raf) & set(maf))
    if len(common) < 2:
        raise ValueError("need paired Figure 5B x_RAF/x_MAF points")
    t = np.asarray(common, dtype=float)
    xr = np.asarray([raf[T] for T in common], dtype=float)
    xm = np.asarray([maf[T] for T in common], dtype=float)
    return t, xr, xm


def build_regularized_oaf_state(
    csv_path: str | Path,
    *,
    temperature_C: float,
    reference_low_C: float = -30.0,
    reference_high_C: float = 40.0,
    crystal_fraction: float = 0.60,
    iaf_fraction: float = 0.20,
    total_oaf_fraction: float = 0.20,
    max_q_disagreement: float = 0.08,
) -> RegularizedOAFDevitrificationState:
    """Return a non-negative OAF subpartition constrained by Figure 5B progress.

    The supported range is intentionally limited to [-30, 40] C, matching the
    temperature window used for the source BDS/OAF dielectric analysis.
    """
    t, xr, xm = _paired_curves(csv_path)
    T = float(temperature_C)
    lo = float(reference_low_C)
    hi = float(reference_high_C)
    if not (lo <= T <= hi):
        raise ValueError(f"temperature must lie in [{lo:g}, {hi:g}] C")
    if lo < t.min() or hi > t.max():
        raise ValueError("reference temperatures are outside digitized Figure 5B")

    xr_T = float(np.interp(T, t, xr))
    xm_T = float(np.interp(T, t, xm))
    xr_lo = float(np.interp(lo, t, xr))
    xm_lo = float(np.interp(lo, t, xm))
    xr_hi = float(np.interp(hi, t, xr))
    xm_hi = float(np.interp(hi, t, xm))

    raf_den = xr_lo - xr_hi
    maf_den = xm_hi - xm_lo
    if raf_den <= 0.0 or maf_den <= 0.0:
        raise ValueError("Figure 5B endpoints do not represent devitrification")

    q_raf = (xr_lo - xr_T) / raf_den
    q_maf = (xm_T - xm_lo) / maf_den
    q_raf = float(np.clip(q_raf, 0.0, 1.0))
    q_maf = float(np.clip(q_maf, 0.0, 1.0))
    disagreement = abs(q_raf - q_maf)
    if disagreement > max_q_disagreement:
        raise ValueError(
            f"RAF-loss and MAF-gain progress disagree by {disagreement:.3f}; "
            "review Figure 5B digitization/source assumptions"
        )
    q = float(0.5 * (q_raf + q_maf))

    fc = float(crystal_fraction)
    fi = float(iaf_fraction)
    fo = float(total_oaf_fraction)
    if min(fc, fi, fo) < 0.0 or not np.isclose(fc + fi + fo, 1.0, atol=1e-12):
        raise ValueError("project regularized structural fractions must sum to one")
    roaf = fo * (1.0 - q)
    moaf = fo * q

    return RegularizedOAFDevitrificationState(
        temperature_C=T,
        x_raf_raw=xr_T,
        x_maf_raw=xm_T,
        q_from_raf_loss=q_raf,
        q_from_maf_gain=q_maf,
        q_devitrification=q,
        q_internal_disagreement=disagreement,
        crystal=fc,
        roaf=float(roaf),
        moaf=float(moaf),
        iaf=fi,
    )
