"""Same-state BOPVDF dielectric-relaxation bridge for v0.1.12.

This module intentionally stops at the response that Rui et al. 2022 can support
without mixing structural states:

* film static permittivity epsilon_c(T) from Figure 3A;
* alpha-relaxation peak frequency f_d(T) from Figure 3B;
* eta_cr=0.52 and epsilon_cr=3 from the article's two-phase treatment;
* the inferred *combined amorphous* static permittivity epsilon_am(T).

The source explicitly states that the BDS/Kirkwood-Frohlich treatment cannot
separate OAF and IAF.  Therefore no OAF-only tau(T), epsilon*(omega,T), or
mobility is created here.

A single-Debye response can be constructed only after the caller supplies an
independent fast-limit amorphous permittivity.  That operation is a numerical
bridge, not a claim that the measured alpha relaxation is exactly Debye.

Complex-permittivity convention: eps* = eps' - i eps'', so dielectric loss is
positive and equals -Im(eps*).
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from pvdf_pf.calibration.bds import debye_complex_permittivity
from pvdf_pf.calibration.source_data import DigitizedPoint, load_digitized_curve_csv

R_GAS = 8.31446261815324  # J mol^-1 K^-1


@dataclass(frozen=True)
class SameStateBDSPoint:
    sample_state: str
    temperature_C: float
    epsilon_film_static: float
    eta_crystal: float
    epsilon_crystal: float
    epsilon_amorphous_static: float
    ln_f_d: float
    f_d_hz: float
    tau_alpha_s: float
    response_assignment: str = "combined_OAF_plus_IAF_amorphous_response"

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class ArrheniusFit:
    sample_state: str
    n_points: int
    min_temperature_C_exclusive: float
    slope_vs_1000_over_T: float
    intercept_ln_hz: float
    activation_energy_kj_mol: float
    r_squared: float

    def to_dict(self) -> dict:
        return asdict(self)


def _curve(
    points: list[DigitizedPoint],
    *,
    sample_state: str,
    quantity: str,
) -> list[DigitizedPoint]:
    sample = str(sample_state).strip().lower()
    selected = [
        p
        for p in points
        if p.sample_state.strip().lower() == sample and p.quantity == quantity
    ]
    selected.sort(key=lambda p: p.temperature_C)
    if len(selected) < 2:
        raise ValueError(
            f"insufficient {quantity!r} points for sample_state={sample!r}"
        )
    return selected


def _exact_value(
    curve: list[DigitizedPoint],
    temperature_C: float,
    *,
    atol: float = 1e-9,
) -> float:
    T = float(temperature_C)
    values = [p.value for p in curve if abs(float(p.temperature_C) - T) <= atol]
    if len(values) != 1:
        raise ValueError(
            f"expected exactly one source point at T={T:g} C, got {len(values)}"
        )
    return float(values[0])


def infer_combined_amorphous_static_permittivity(
    epsilon_film_static: float,
    *,
    eta_crystal: float = 0.52,
    epsilon_crystal: float = 3.0,
) -> float:
    r"""Apply the article's parallel two-phase closure.

    epsilon_film = eta_cr * epsilon_cr + (1-eta_cr) * epsilon_am.

    `epsilon_am` is the combined amorphous response.  It is not an OAF-only
    permittivity.
    """
    eps_film = float(epsilon_film_static)
    eta = float(eta_crystal)
    eps_cr = float(epsilon_crystal)
    if eps_film <= 0.0 or eps_cr <= 0.0:
        raise ValueError("permittivities must be positive")
    if not (0.0 <= eta < 1.0):
        raise ValueError("eta_crystal must satisfy 0 <= eta_crystal < 1")
    eps_am = (eps_film - eta * eps_cr) / (1.0 - eta)
    if eps_am <= 0.0:
        raise ValueError(
            "two-phase inversion produced non-positive combined amorphous permittivity"
        )
    return float(eps_am)


def build_same_state_bds_point(
    rui2022_dir: str | Path,
    *,
    sample_state: str,
    temperature_C: float,
    eta_crystal: float = 0.52,
    epsilon_crystal: float = 3.0,
) -> SameStateBDSPoint:
    root = Path(rui2022_dir)
    sample = str(sample_state).strip().lower()
    if sample not in {"unpoled", "poled"}:
        raise ValueError("sample_state must be 'unpoled' or 'poled'")

    film_points = load_digitized_curve_csv(root / "figure_3a_film_eps_digitized.csv")
    fd_points = load_digitized_curve_csv(root / "figure_3b_lnfd_digitized.csv")

    eps_film = _exact_value(
        _curve(film_points, sample_state=sample, quantity="epsilon_c_film"),
        temperature_C,
    )
    ln_fd = _exact_value(
        _curve(fd_points, sample_state=sample, quantity="ln_f_d"),
        temperature_C,
    )
    fd = float(np.exp(ln_fd))
    tau = float(1.0 / (2.0 * np.pi * fd))
    eps_am = infer_combined_amorphous_static_permittivity(
        eps_film,
        eta_crystal=eta_crystal,
        epsilon_crystal=epsilon_crystal,
    )
    return SameStateBDSPoint(
        sample_state=sample,
        temperature_C=float(temperature_C),
        epsilon_film_static=float(eps_film),
        eta_crystal=float(eta_crystal),
        epsilon_crystal=float(epsilon_crystal),
        epsilon_amorphous_static=eps_am,
        ln_f_d=float(ln_fd),
        f_d_hz=fd,
        tau_alpha_s=tau,
    )


def same_state_bds_series(
    rui2022_dir: str | Path,
    *,
    sample_state: str,
    temperatures_C: tuple[float, ...] = (
        -30.0,
        -20.0,
        -10.0,
        0.0,
        10.0,
        20.0,
        30.0,
        40.0,
    ),
    eta_crystal: float = 0.52,
    epsilon_crystal: float = 3.0,
) -> list[SameStateBDSPoint]:
    return [
        build_same_state_bds_point(
            rui2022_dir,
            sample_state=sample_state,
            temperature_C=T,
            eta_crystal=eta_crystal,
            epsilon_crystal=epsilon_crystal,
        )
        for T in temperatures_C
    ]


def fit_high_temperature_arrhenius(
    points: list[SameStateBDSPoint],
    *,
    min_temperature_C_exclusive: float = 0.0,
) -> ArrheniusFit:
    """Fit ln(f_d) = intercept + slope*(1000/T) above the regime boundary."""
    selected = [
        p for p in points if p.temperature_C > float(min_temperature_C_exclusive)
    ]
    if len(selected) < 3:
        raise ValueError("at least three high-temperature points are required")
    sample_states = {p.sample_state for p in selected}
    if len(sample_states) != 1:
        raise ValueError("Arrhenius fit requires one sample_state")

    temperature_K = np.asarray([p.temperature_C + 273.15 for p in selected])
    x = 1000.0 / temperature_K
    y = np.asarray([p.ln_f_d for p in selected])
    slope, intercept = np.polyfit(x, y, 1)
    yhat = slope * x + intercept
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 if ss_tot == 0.0 else 1.0 - ss_res / ss_tot

    # Because x=1000/T, slope=-Ea/(1000 R).  Converting Ea to kJ/mol gives
    # Ea[kJ/mol] = -slope * R[J/mol/K].
    ea_kj_mol = float(-slope * R_GAS)
    return ArrheniusFit(
        sample_state=next(iter(sample_states)),
        n_points=len(selected),
        min_temperature_C_exclusive=float(min_temperature_C_exclusive),
        slope_vs_1000_over_T=float(slope),
        intercept_ln_hz=float(intercept),
        activation_energy_kj_mol=ea_kj_mol,
        r_squared=float(r2),
    )


def amorphous_complex_permittivity_debye_bridge(
    frequency_hz: np.ndarray | float,
    point: SameStateBDSPoint,
    *,
    epsilon_amorphous_fast: float,
) -> np.ndarray:
    """Construct a one-relaxation complex epsilon_am* after supplying eps_am,fast.

    `epsilon_amorphous_fast` is intentionally not inferred in v0.1.12.  It must
    come from an independently digitized high-frequency plateau or another
    justified measurement/model.  The returned spectrum is therefore a bridge
    for numerical coupling, not a direct OAF/IAF-resolved fit.
    """
    eps_fast = float(epsilon_amorphous_fast)
    eps_static = float(point.epsilon_amorphous_static)
    if not (0.0 < eps_fast <= eps_static):
        raise ValueError(
            "require 0 < epsilon_amorphous_fast <= epsilon_amorphous_static"
        )
    return debye_complex_permittivity(
        frequency_hz,
        eps_inf=eps_fast,
        delta_eps=eps_static - eps_fast,
        tau_s=point.tau_alpha_s,
    )


def film_complex_permittivity_parallel_bridge(
    frequency_hz: np.ndarray | float,
    point: SameStateBDSPoint,
    *,
    epsilon_amorphous_fast: float,
) -> np.ndarray:
    """Parallel crystal + combined-amorphous spectrum consistent with Eq. 20."""
    eps_am = amorphous_complex_permittivity_debye_bridge(
        frequency_hz,
        point,
        epsilon_amorphous_fast=epsilon_amorphous_fast,
    )
    return (
        point.eta_crystal * point.epsilon_crystal
        + (1.0 - point.eta_crystal) * eps_am
    )
