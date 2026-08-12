"""Audit the sample-state transfer behind Rui-2022 Figure S2.

The 2022 paper contains two distinct material states:

* the BDS target is unpoled/poled biaxially oriented PVDF (BOPVDF), for which
  the main text states crystallinity x_c ~= 0.52;
* the RAF/MAF mobility donor in Figure 5 is a melt-recrystallized BOPVDF sample
  with alpha crystals and x_c = 0.59 (rounded to eta_cr=0.60 in SI Sec. S2).

Rui-2021 independently reports, for highly poled BOPVDF, x_c ~= 0.52 and an
OAF content ~= 0.28.  The remaining ~=0.20 is therefore used here only as a
closed same-state *structural fraction proxy*.  These approximate source
fractions are not promoted as exact local volume fractions.

This module asks a limited inverse question: under the same-state closed proxy
(0.52 crystal + 0.28 OAF + 0.20 IAF), what effective OAF small-signal response
would be required to reproduce the measured BOPVDF film permittivity?

The answer is a film-closure effective response, not a measured intrinsic OAF
permittivity and not a TDGL background dielectric constant.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from pvdf_pf.calibration.iaf_permittivity import epsilon_iaf_from_melt_extrapolation
from pvdf_pf.calibration.source_data import DigitizedPoint, load_digitized_curve_csv


@dataclass(frozen=True)
class SourceStateDefinition:
    name: str
    sample_state: str
    crystal_fraction: float
    oaf_fraction: float | None
    iaf_fraction: float | None
    raf_maf_total: float | None
    role: str
    provenance: str

    @property
    def closed_sum(self) -> float | None:
        if self.oaf_fraction is not None and self.iaf_fraction is not None:
            return float(self.crystal_fraction + self.oaf_fraction + self.iaf_fraction)
        if self.raf_maf_total is not None:
            return float(self.crystal_fraction + self.raf_maf_total)
        return None

    def to_dict(self) -> dict:
        out = asdict(self)
        out["closed_sum"] = self.closed_sum
        return out


@dataclass(frozen=True)
class ClosedOAFResponsePoint:
    sample_state: str
    temperature_C: float
    epsilon_film_measured: float
    epsilon_crystal: float
    epsilon_iaf_kf: float
    fraction_crystal: float
    fraction_oaf_proxy: float
    fraction_iaf_proxy: float
    epsilon_oaf_effective_required: float
    epsilon_moaf_figure_s2: float
    effective_oaf_minus_figure_s2_moaf: float
    provenance: str = "SAME_STATE_CLOSED_STRUCTURAL_PROXY_INVERSION"

    def to_dict(self) -> dict:
        return asdict(self)


def source_state_definitions() -> tuple[SourceStateDefinition, ...]:
    return (
        SourceStateDefinition(
            name="rui2022_bds_target",
            sample_state="unpoled/poled BOPVDF",
            crystal_fraction=0.52,
            oaf_fraction=None,
            iaf_fraction=None,
            raf_maf_total=0.48,
            role="film dielectric target state",
            provenance="Rui-2022 main text: both unpoled and poled BOPVDF crystallinity around 0.52",
        ),
        SourceStateDefinition(
            name="rui2022_tmdsc_mobility_donor",
            sample_state="melt-recrystallized BOPVDF, alpha phase",
            crystal_fraction=0.59,
            oaf_fraction=None,
            iaf_fraction=None,
            raf_maf_total=0.41,
            role="RAF/MAF devitrification donor state",
            provenance="Rui-2022 Figure 5: x_c=0.59, x_RAF+x_MAF approximately 0.41",
        ),
        SourceStateDefinition(
            name="rui2022_si_s2_rounded_model",
            sample_state="SI transfer model",
            crystal_fraction=0.60,
            oaf_fraction=0.20,
            iaf_fraction=0.20,
            raf_maf_total=None,
            role="Figure S2 local-MOAF inverse model",
            provenance="Rui-2022 SI S2: eta_cr=0.60, eta_IAF approximately eta_OAF=0.20",
        ),
        SourceStateDefinition(
            name="rui2021_poled_bopvdf_structural_proxy",
            sample_state="highly poled BOPVDF",
            crystal_fraction=0.52,
            oaf_fraction=0.28,
            iaf_fraction=0.20,
            raf_maf_total=None,
            role="same-state closed structural proxy for film inversion",
            provenance=(
                "Rui-2021 reports crystallinity approximately 0.52 and OAF content approximately 0.28; "
                "IAF=0.20 is closure remainder, not an independent local-volume measurement"
            ),
        ),
    )


def _curve(points: list[DigitizedPoint], *, sample_state: str, quantity: str) -> list[DigitizedPoint]:
    selected = [
        p for p in points
        if p.sample_state.lower() == sample_state.lower() and p.quantity == quantity
    ]
    selected.sort(key=lambda p: p.temperature_C)
    if len(selected) < 2:
        raise ValueError(f"insufficient {quantity} points for {sample_state!r}")
    return selected


def _exact(curve: list[DigitizedPoint], temperature_C: float, atol: float = 1e-9) -> float:
    values = [p.value for p in curve if abs(p.temperature_C - float(temperature_C)) <= atol]
    if len(values) != 1:
        raise ValueError(f"expected one source value at {temperature_C:g} C, got {len(values)}")
    return float(values[0])


def closed_oaf_effective_response(
    rui2022_dir: str | Path,
    *,
    sample_state: str,
    temperature_C: float,
    fraction_crystal: float = 0.52,
    fraction_oaf: float = 0.28,
    fraction_iaf: float = 0.20,
    epsilon_crystal: float = 3.0,
) -> ClosedOAFResponsePoint:
    """Invert an effective OAF response using a closed same-BOPVDF fraction proxy."""
    total = float(fraction_crystal + fraction_oaf + fraction_iaf)
    if not np.isclose(total, 1.0, atol=1e-12):
        raise ValueError(f"closed structural proxy must sum to one, got {total}")
    if min(fraction_crystal, fraction_oaf, fraction_iaf) < 0.0:
        raise ValueError("fractions must be non-negative")
    if fraction_oaf <= 0.0:
        raise ValueError("OAF fraction must be positive for inversion")

    root = Path(rui2022_dir)
    sample = str(sample_state).strip().lower()
    T = float(temperature_C)

    film = load_digitized_curve_csv(root / "figure_3a_film_eps_digitized.csv")
    eps_film = _exact(_curve(film, sample_state=sample, quantity="epsilon_c_film"), T)
    eps_iaf, _ = epsilon_iaf_from_melt_extrapolation(
        root / "figure_1b_melt_eps_digitized.csv",
        temperature_C=T,
    )

    s2 = load_digitized_curve_csv(root / "figure_s2_digitized.csv")
    eps_moaf_s2 = _exact(_curve(s2, sample_state=sample, quantity="epsilon_MOAF"), T)

    eps_oaf_eff = float(
        (
            eps_film
            - fraction_crystal * float(epsilon_crystal)
            - fraction_iaf * eps_iaf
        )
        / fraction_oaf
    )

    return ClosedOAFResponsePoint(
        sample_state=sample,
        temperature_C=T,
        epsilon_film_measured=eps_film,
        epsilon_crystal=float(epsilon_crystal),
        epsilon_iaf_kf=float(eps_iaf),
        fraction_crystal=float(fraction_crystal),
        fraction_oaf_proxy=float(fraction_oaf),
        fraction_iaf_proxy=float(fraction_iaf),
        epsilon_oaf_effective_required=eps_oaf_eff,
        epsilon_moaf_figure_s2=float(eps_moaf_s2),
        effective_oaf_minus_figure_s2_moaf=float(eps_oaf_eff - eps_moaf_s2),
    )


def closed_oaf_series(
    rui2022_dir: str | Path,
    *,
    sample_state: str,
    temperatures_C: tuple[float, ...] = (-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0),
) -> list[ClosedOAFResponsePoint]:
    return [
        closed_oaf_effective_response(
            rui2022_dir,
            sample_state=sample_state,
            temperature_C=T,
        )
        for T in temperatures_C
    ]


def summarize_closed_oaf(points: list[ClosedOAFResponsePoint]) -> dict:
    if not points:
        raise ValueError("points must not be empty")
    eps = np.asarray([p.epsilon_oaf_effective_required for p in points], dtype=float)
    s2 = np.asarray([p.epsilon_moaf_figure_s2 for p in points], dtype=float)
    T = np.asarray([p.temperature_C for p in points], dtype=float)
    slope, intercept = np.polyfit(T, eps, 1)
    return {
        "n_points": len(points),
        "epsilon_oaf_effective_min": float(eps.min()),
        "epsilon_oaf_effective_max": float(eps.max()),
        "epsilon_oaf_effective_at_lowest_T": float(eps[np.argmin(T)]),
        "epsilon_oaf_effective_at_highest_T": float(eps[np.argmax(T)]),
        "linear_slope_per_C_diagnostic": float(slope),
        "linear_intercept_diagnostic": float(intercept),
        "rmse_vs_figure_s2_moaf": float(np.sqrt(np.mean((eps - s2) ** 2))),
        "interpretation": (
            "epsilon_OAF,eff is the response required by same-state film closure. "
            "Figure-S2 epsilon_MOAF is a different local/mobile-OAF inverse quantity and need not coincide."
        ),
    }
