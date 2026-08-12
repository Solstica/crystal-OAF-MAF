"""IAF permittivity extrapolation from Rui et al. 2022 Figure 1B.

Supporting Information Section S2 states that epsilon_IAF(T) is obtained by
extrapolating the molten-PVDF permittivity above about 147 C to low temperature.
This module implements exactly that source-model step from the traceably digitized
main-text Figure 1B data.

Important scope
---------------
The resulting epsilon_IAF(T) is a SOURCE_MODEL_EXTRAPOLATION used for the
small-signal dielectric reconstruction. It is not a direct local measurement and
must not be copied into the TDGL background permittivity `eps_b` without a
separate double-counting analysis.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from pvdf_pf.calibration.source_data import load_digitized_curve_csv


@dataclass(frozen=True)
class MeltPermittivityLinearFit:
    slope_per_C: float
    intercept: float
    r_squared: float
    rmse: float
    temperature_min_C: float
    temperature_max_C: float
    n_points: int
    provenance: str = "SOURCE_MODEL_EXTRAPOLATION"

    def epsilon(self, temperature_C: float) -> float:
        value = self.slope_per_C * float(temperature_C) + self.intercept
        if value <= 0.0:
            raise ValueError("linear melt extrapolation produced non-positive permittivity")
        return float(value)

    def to_dict(self) -> dict:
        return asdict(self)


def fit_molten_pvdf_permittivity(csv_path: str | Path) -> MeltPermittivityLinearFit:
    """Fit epsilon_c(T)=a*T+b to the digitized molten-PVDF Figure 1B data."""
    points = load_digitized_curve_csv(csv_path)
    selected = [
        p for p in points
        if p.figure == "Figure1"
        and p.panel == "B"
        and p.sample_state.lower() == "molten_pvdf"
        and p.quantity == "epsilon_c_melt"
    ]
    selected.sort(key=lambda p: p.temperature_C)
    if len(selected) < 3:
        raise ValueError("need at least three Figure 1B molten-PVDF points")

    T = np.asarray([p.temperature_C for p in selected], dtype=float)
    eps = np.asarray([p.value for p in selected], dtype=float)
    if np.unique(T).size != T.size:
        raise ValueError("duplicate Figure 1B temperatures")
    if np.any(eps <= 0.0):
        raise ValueError("permittivity must be positive")

    slope, intercept = np.polyfit(T, eps, 1)
    pred = slope * T + intercept
    residual = eps - pred
    ss_res = float(np.sum(residual**2))
    ss_tot = float(np.sum((eps - np.mean(eps))**2))
    r_squared = 1.0 if ss_tot == 0.0 else 1.0 - ss_res / ss_tot
    rmse = float(np.sqrt(np.mean(residual**2)))

    return MeltPermittivityLinearFit(
        slope_per_C=float(slope),
        intercept=float(intercept),
        r_squared=float(r_squared),
        rmse=rmse,
        temperature_min_C=float(T.min()),
        temperature_max_C=float(T.max()),
        n_points=int(T.size),
    )


def epsilon_iaf_from_melt_extrapolation(
    csv_path: str | Path,
    *,
    temperature_C: float,
    allowed_low_C: float = -30.0,
    allowed_high_C: float = 40.0,
) -> tuple[float, MeltPermittivityLinearFit]:
    """Return source-model epsilon_IAF(T) in the project BDS window.

    The SI explicitly requests extrapolation from the melt. The project restricts
    use to -30..40 C because that is the temperature window of the extracted
    MOAF/BDS source states; it does not claim the linear law is universal.
    """
    T = float(temperature_C)
    if not (allowed_low_C <= T <= allowed_high_C):
        raise ValueError(
            f"epsilon_IAF source-model extrapolation is restricted to "
            f"[{allowed_low_C:g}, {allowed_high_C:g}] C"
        )
    fit = fit_molten_pvdf_permittivity(csv_path)
    return fit.epsilon(T), fit
