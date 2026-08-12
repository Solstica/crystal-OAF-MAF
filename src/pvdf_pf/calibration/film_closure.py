"""v0.1.8 closure audit against Rui et al. 2022 main-text Figure 3A.

The purpose of this module is diagnostic. It asks whether the local dielectric
quantities currently carried by the project reproduce the measured BOPVDF film
small-signal permittivity under the same ideal parallel-capacitor accounting used
by the source paper.

No parameter is refitted here. A closure failure is retained as information rather
than hidden by changing epsilon_MOAF, epsilon_IAF, or the project devitrification
mapping.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from pvdf_pf.calibration.four_phase_dielectric import build_four_phase_permittivity_state
from pvdf_pf.calibration.project_devitrification import build_regularized_oaf_state
from pvdf_pf.calibration.source_data import DigitizedPoint, load_digitized_curve_csv


@dataclass(frozen=True)
class FilmClosurePoint:
    sample_state: str
    temperature_C: float
    epsilon_film_measured: float
    epsilon_film_project_parallel: float
    residual_project: float
    relative_error_project: float
    epsilon_oaf_mean_project: float
    epsilon_oaf_mean_required_for_closure: float
    oaf_mean_gap: float
    eta_roaf_project: float
    eta_moaf_project: float
    epsilon_roaf: float
    epsilon_moaf: float
    epsilon_iaf: float
    literal_si_available: bool
    literal_si_fraction_sum: float | None
    literal_si_eta_roaf: float | None
    literal_si_eta_moaf: float | None
    epsilon_film_literal_si: float | None
    residual_literal_si: float | None

    def to_dict(self) -> dict:
        return asdict(self)


def _curve(points: list[DigitizedPoint], *, sample_state: str, quantity: str) -> list[DigitizedPoint]:
    sample = str(sample_state).strip().lower()
    selected = [
        p for p in points
        if p.sample_state.lower() == sample and p.quantity == quantity
    ]
    selected.sort(key=lambda p: p.temperature_C)
    if len(selected) < 2:
        raise ValueError(f"insufficient {quantity} points for sample_state={sample!r}")
    return selected


def _exact_point(curve: list[DigitizedPoint], temperature_C: float, *, atol: float = 1e-9) -> float:
    T = float(temperature_C)
    matches = [p.value for p in curve if abs(float(p.temperature_C) - T) <= atol]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one source point at T={T:g} C, got {len(matches)}")
    return float(matches[0])


def audit_film_closure(
    rui2022_dir: str | Path,
    *,
    sample_state: str,
    temperature_C: float,
    eta_crystal: float = 0.60,
    eta_iaf: float = 0.20,
    eta_total_oaf: float = 0.20,
) -> FilmClosurePoint:
    """Evaluate project and literal-SI closure at one source temperature.

    `literal_si_*` implements the algebra written in SI Eq. S4 only when the
    implied eta_MOAF = x_MAF - eta_IAF is non-negative. It is reported as an
    audit, not promoted into the project state; the raw Figure 5B closure is about
    0.41 because the main-text calorimetry used x_c=0.59.
    """
    root = Path(rui2022_dir)
    sample = str(sample_state).strip().lower()
    T = float(temperature_C)

    film_points = load_digitized_curve_csv(root / "figure_3a_film_eps_digitized.csv")
    eps_meas = _exact_point(
        _curve(film_points, sample_state=sample, quantity="epsilon_c_film"), T
    )

    phase_state = build_four_phase_permittivity_state(
        root,
        temperature_C=T,
        sample_state=sample,
    )
    dev = build_regularized_oaf_state(root / "figure_5b_maintext_digitized.csv", temperature_C=T)
    if not np.isclose(dev.crystal, eta_crystal, atol=1e-12):
        raise ValueError("closure audit eta_crystal differs from v0.1.6 structural state")
    if not np.isclose(dev.iaf, eta_iaf, atol=1e-12):
        raise ValueError("closure audit eta_iaf differs from v0.1.6 structural state")
    if not np.isclose(dev.total_oaf, eta_total_oaf, atol=1e-12):
        raise ValueError("closure audit eta_total_oaf differs from v0.1.6 structural state")
    eta_roaf = float(dev.roaf)
    eta_moaf = float(dev.moaf)

    eps_project = float(
        eta_crystal * phase_state.epsilon_crystal
        + eta_roaf * phase_state.epsilon_roaf
        + eta_moaf * phase_state.epsilon_moaf
        + eta_iaf * phase_state.epsilon_iaf
    )
    eps_oaf_project = float(
        (eta_roaf * phase_state.epsilon_roaf + eta_moaf * phase_state.epsilon_moaf)
        / eta_total_oaf
    )
    eps_oaf_required = float(
        (
            eps_meas
            - eta_crystal * phase_state.epsilon_crystal
            - eta_iaf * phase_state.epsilon_iaf
        )
        / eta_total_oaf
    )

    fig5 = load_digitized_curve_csv(root / "figure_5b_maintext_digitized.csv")
    x_raf = _exact_point(
        _curve(fig5, sample_state="melt-recrystallized", quantity="x_RAF"), T
    )
    x_maf = _exact_point(
        _curve(fig5, sample_state="melt-recrystallized", quantity="x_MAF"), T
    )
    literal_eta_roaf = float(x_raf)
    literal_eta_moaf = float(x_maf - eta_iaf)
    literal_sum = float(eta_crystal + eta_iaf + literal_eta_roaf + literal_eta_moaf)
    literal_available = literal_eta_moaf >= -1e-12

    eps_literal = None
    residual_literal = None
    if literal_available:
        literal_eta_moaf = max(literal_eta_moaf, 0.0)
        eps_literal = float(
            eta_crystal * phase_state.epsilon_crystal
            + literal_eta_roaf * phase_state.epsilon_roaf
            + literal_eta_moaf * phase_state.epsilon_moaf
            + eta_iaf * phase_state.epsilon_iaf
        )
        residual_literal = float(eps_literal - eps_meas)

    return FilmClosurePoint(
        sample_state=sample,
        temperature_C=T,
        epsilon_film_measured=eps_meas,
        epsilon_film_project_parallel=eps_project,
        residual_project=float(eps_project - eps_meas),
        relative_error_project=float((eps_project - eps_meas) / eps_meas),
        epsilon_oaf_mean_project=eps_oaf_project,
        epsilon_oaf_mean_required_for_closure=eps_oaf_required,
        oaf_mean_gap=float(eps_oaf_required - eps_oaf_project),
        eta_roaf_project=eta_roaf,
        eta_moaf_project=eta_moaf,
        epsilon_roaf=float(phase_state.epsilon_roaf),
        epsilon_moaf=float(phase_state.epsilon_moaf),
        epsilon_iaf=float(phase_state.epsilon_iaf),
        literal_si_available=bool(literal_available),
        literal_si_fraction_sum=(literal_sum if literal_available else None),
        literal_si_eta_roaf=(literal_eta_roaf if literal_available else None),
        literal_si_eta_moaf=(literal_eta_moaf if literal_available else None),
        epsilon_film_literal_si=eps_literal,
        residual_literal_si=residual_literal,
    )


def closure_series(
    rui2022_dir: str | Path,
    *,
    sample_state: str,
    temperatures_C: tuple[float, ...] = (-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0),
) -> list[FilmClosurePoint]:
    return [
        audit_film_closure(
            rui2022_dir,
            sample_state=sample_state,
            temperature_C=T,
        )
        for T in temperatures_C
    ]


def summarize_closure(points: list[FilmClosurePoint]) -> dict[str, float | int]:
    if not points:
        raise ValueError("points must not be empty")
    residual = np.asarray([p.residual_project for p in points], dtype=float)
    rel = np.asarray([p.relative_error_project for p in points], dtype=float)
    gaps = np.asarray([p.oaf_mean_gap for p in points], dtype=float)
    return {
        "n_points": int(len(points)),
        "rmse_project": float(np.sqrt(np.mean(residual**2))),
        "mae_project": float(np.mean(np.abs(residual))),
        "mean_relative_error_project": float(np.mean(rel)),
        "max_abs_relative_error_project": float(np.max(np.abs(rel))),
        "mean_required_oaf_gap": float(np.mean(gaps)),
        "max_required_oaf_gap": float(np.max(gaps)),
    }
