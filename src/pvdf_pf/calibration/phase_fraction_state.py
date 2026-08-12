"""Reconstruct crystal/ROAF/MOAF/IAF fractions from Rui et al. 2022.

The Supporting Information defines, for the source-model accounting,

    eta_cr + eta_ROAF + eta_MOAF + eta_IAF = 1

and uses eta_cr = 0.6 and eta_IAF = 0.2.  For temperatures above Tg,

    x_RAF = eta_ROAF
    x_MAF = eta_MOAF + eta_IAF.

This module encodes that algebra only.  It does not infer x_RAF(T) or x_MAF(T);
those curves must be digitized from Figure 5B of the main article.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from pvdf_pf.calibration.source_data import DigitizedPoint, load_digitized_curve_csv
from pvdf_pf.calibration.oaf_partition import FourFractionState


@dataclass(frozen=True)
class SourcePhaseFractionState:
    sample_state: str
    temperature_C: float
    x_raf: float
    x_maf: float
    crystal: float
    roaf: float
    moaf: float
    iaf: float
    mobile_fraction_within_oaf: float

    def to_dict(self) -> dict:
        return asdict(self)

    @property
    def total_oaf(self) -> float:
        return float(self.roaf + self.moaf)

    def as_four_fraction_state(self) -> FourFractionState:
        state = FourFractionState(
            crystal=self.crystal,
            roaf=self.roaf,
            moaf=self.moaf,
            iaf=self.iaf,
        )
        state.validate()
        return state


def reconstruct_above_tg(
    *,
    x_raf: float,
    x_maf: float,
    temperature_C: float,
    sample_state: str,
    crystal_fraction: float = 0.6,
    iaf_fraction: float = 0.2,
    atol: float = 5e-3,
) -> SourcePhaseFractionState:
    """Apply SI Eq. S4 to one above-Tg RAF/MAF pair.

    `x_raf` and `x_maf` are DSC/TMDSC phase fractions from the main-text
    Figure 5B.  The reconstruction is source-model-specific and should not be
    mixed with the separate 0.52/0.28/0.20 three-phase parameter set.
    """
    xr = float(x_raf)
    xm = float(x_maf)
    fc = float(crystal_fraction)
    fi = float(iaf_fraction)
    vals = np.asarray([xr, xm, fc, fi, temperature_C], dtype=float)
    if np.any(~np.isfinite(vals)):
        raise ValueError("all fractions and temperature must be finite")
    if min(xr, xm, fc, fi) < 0.0:
        raise ValueError("fractions must be non-negative")
    if fc >= 1.0 or fi >= 1.0:
        raise ValueError("crystal and IAF fractions must be below one")

    amorphous_total = 1.0 - fc
    if not np.isclose(xr + xm, amorphous_total, atol=atol):
        raise ValueError(
            "source RAF + MAF must close to 1 - crystal; "
            f"got {xr + xm:.6g} vs {amorphous_total:.6g}"
        )
    if xm < fi - atol:
        raise ValueError(
            "above Tg, SI Eq. S4 requires x_MAF = eta_MOAF + eta_IAF; "
            "x_MAF cannot be smaller than the assumed IAF fraction"
        )

    roaf = xr
    moaf = max(0.0, xm - fi)
    total_oaf = roaf + moaf
    expected_oaf = 1.0 - fc - fi
    if not np.isclose(total_oaf, expected_oaf, atol=atol):
        raise ValueError(
            "reconstructed ROAF + MOAF does not match total OAF implied by "
            "crystal and IAF fractions"
        )
    mobile_within_oaf = 0.0 if total_oaf <= atol else moaf / total_oaf

    state = SourcePhaseFractionState(
        sample_state=str(sample_state).strip().lower(),
        temperature_C=float(temperature_C),
        x_raf=xr,
        x_maf=xm,
        crystal=fc,
        roaf=roaf,
        moaf=moaf,
        iaf=fi,
        mobile_fraction_within_oaf=float(mobile_within_oaf),
    )
    state.as_four_fraction_state()
    return state


def below_tg_limit(
    *,
    temperature_C: float,
    sample_state: str,
    crystal_fraction: float = 0.6,
    iaf_fraction: float = 0.2,
) -> SourcePhaseFractionState:
    """Return the SI limiting partition before OAF devitrification.

    Below Tg the source states x_MAF = 0 and assigns the whole noncrystalline
    material to RAF for DSC accounting.  Spatially, the present project keeps
    IAF as a separate structural region and places the entire OAF in ROAF.
    """
    fc = float(crystal_fraction)
    fi = float(iaf_fraction)
    fo = 1.0 - fc - fi
    if min(fc, fi, fo) < 0.0:
        raise ValueError("invalid fixed source fractions")
    state = SourcePhaseFractionState(
        sample_state=str(sample_state).strip().lower(),
        temperature_C=float(temperature_C),
        x_raf=float(fo + fi),
        x_maf=0.0,
        crystal=fc,
        roaf=fo,
        moaf=0.0,
        iaf=fi,
        mobile_fraction_within_oaf=0.0,
    )
    state.as_four_fraction_state()
    return state


def _select_curve(
    points: list[DigitizedPoint], *, sample_state: str, quantity: str
) -> list[DigitizedPoint]:
    sample = str(sample_state).strip().lower()
    selected = [
        p for p in points if p.sample_state.lower() == sample and p.quantity == quantity
    ]
    selected.sort(key=lambda p: p.temperature_C)
    if len(selected) < 2:
        raise ValueError(
            f"need at least two Figure 5B points for {sample!r}, quantity={quantity!r}"
        )
    return selected


def reconstruct_curve_from_figure5b(
    csv_path: str | Path,
    *,
    sample_state: str,
    crystal_fraction: float = 0.6,
    iaf_fraction: float = 0.2,
    atol: float = 5e-3,
) -> list[SourcePhaseFractionState]:
    """Reconstruct all temperatures common to digitized x_RAF/x_MAF curves."""
    points = load_digitized_curve_csv(csv_path)
    raf = _select_curve(points, sample_state=sample_state, quantity="x_RAF")
    maf = _select_curve(points, sample_state=sample_state, quantity="x_MAF")
    raf_by_t = {float(p.temperature_C): float(p.value) for p in raf}
    maf_by_t = {float(p.temperature_C): float(p.value) for p in maf}
    common = sorted(set(raf_by_t) & set(maf_by_t))
    if len(common) < 2:
        raise ValueError("x_RAF and x_MAF need at least two common temperatures")
    return [
        reconstruct_above_tg(
            x_raf=raf_by_t[T],
            x_maf=maf_by_t[T],
            temperature_C=T,
            sample_state=sample_state,
            crystal_fraction=crystal_fraction,
            iaf_fraction=iaf_fraction,
            atol=atol,
        )
        for T in common
    ]
