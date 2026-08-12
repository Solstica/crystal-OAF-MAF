"""Rui-2021/2022 amorphous-PVDF dielectric extrapolation for the IAF.

The 2022 Supporting Information says that epsilon_IAF(T) is obtained by
extrapolation from molten PVDF and explicitly points to Figure 9A of Rui et al.,
J. Mater. Chem. C 2021, 9, 894-907.  The 2021 SI does *not* linearly extrapolate
permittivity itself.  It:

1. converts measured molten-PVDF epsilon_s(T) to the Kirkwood-Frohlich g(T),
2. fits g(T) linearly versus absolute temperature,
3. extrapolates g(T) to lower temperature, and
4. solves the Kirkwood-Frohlich equation for epsilon_s(T).

This module implements that source-model construction using the traceably
digitized molten-PVDF Figure 1B data already carried by the project.

Important scope
---------------
The result is a source-model estimate of amorphous/IAF small-signal
permittivity.  It is not a direct local measurement and must not be copied into
the TDGL background permittivity ``eps_b`` without a separate double-counting
analysis.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from pvdf_pf.calibration.source_data import load_digitized_curve_csv


# Exact SI/input values used by Rui 2021 Sec. IV unless noted otherwise.
EPSILON_0_F_PER_M = 8.8541878128e-12
BOLTZMANN_J_PER_K = 1.380649e-23
AVOGADRO_PER_MOL = 6.02214076e23
DEBYE_C_M = 3.33564e-30


@dataclass(frozen=True)
class KirkwoodFrohlichMeltFit:
    epsilon_inf: float
    density_kg_per_m3: float
    dipole_moment_D: float
    molar_mass_kg_per_mol: float
    g_slope_per_K: float
    g_intercept: float
    g_r_squared: float
    g_rmse: float
    temperature_min_C: float
    temperature_max_C: float
    n_points: int
    provenance: str = "RUI2021_KIRKWOOD_FROHLICH_SOURCE_MODEL_EXTRAPOLATION"

    def g(self, temperature_C: float) -> float:
        T_K = float(temperature_C) + 273.15
        if T_K <= 0.0:
            raise ValueError("absolute temperature must be positive")
        return float(self.g_slope_per_K * T_K + self.g_intercept)

    def epsilon(self, temperature_C: float) -> float:
        """Solve the positive Kirkwood-Frohlich root for epsilon_s(T)."""
        T_K = float(temperature_C) + 273.15
        if T_K <= 0.0:
            raise ValueError("absolute temperature must be positive")
        g = self.g(temperature_C)
        if g <= 0.0:
            raise ValueError("extrapolated Kirkwood-Frohlich g-factor is non-positive")

        mu = self.dipole_moment_D * DEBYE_C_M
        rhs = (
            AVOGADRO_PER_MOL
            * self.density_kg_per_m3
            * g
            * mu**2
            / (
                9.0
                * EPSILON_0_F_PER_M
                * self.molar_mass_kg_per_mol
                * BOLTZMANN_J_PER_K
                * T_K
            )
        )
        # (eps-eps_inf)(2 eps+eps_inf) / [eps (eps_inf+2)^2] = rhs
        # -> 2 eps^2 - B eps - eps_inf^2 = 0
        eps_inf = self.epsilon_inf
        B = eps_inf + rhs * (eps_inf + 2.0) ** 2
        eps = (B + np.sqrt(B * B + 8.0 * eps_inf * eps_inf)) / 4.0
        if eps <= 0.0:
            raise ValueError("Kirkwood-Frohlich inversion produced non-positive permittivity")
        return float(eps)

    def to_dict(self) -> dict:
        out = asdict(self)
        # Published Rui-2021 anchors are useful executable checks of the source transfer.
        out["source_anchor_epsilon_25C"] = self.epsilon(25.0)
        out["source_anchor_epsilon_127C"] = self.epsilon(127.0)
        return out


def _kf_reduced_lhs(epsilon_s: float, epsilon_inf: float) -> float:
    eps = float(epsilon_s)
    eps_inf = float(epsilon_inf)
    if eps <= 0.0 or eps_inf <= 0.0:
        raise ValueError("permittivities must be positive")
    return float(
        (eps - eps_inf) * (2.0 * eps + eps_inf)
        / (eps * (eps_inf + 2.0) ** 2)
    )


def _g_from_epsilon(
    epsilon_s: float,
    temperature_C: float,
    *,
    epsilon_inf: float,
    density_kg_per_m3: float,
    dipole_moment_D: float,
    molar_mass_kg_per_mol: float,
) -> float:
    T_K = float(temperature_C) + 273.15
    if T_K <= 0.0:
        raise ValueError("absolute temperature must be positive")
    mu = float(dipole_moment_D) * DEBYE_C_M
    reduced = _kf_reduced_lhs(epsilon_s, epsilon_inf)
    return float(
        reduced
        * 9.0
        * EPSILON_0_F_PER_M
        * float(molar_mass_kg_per_mol)
        * BOLTZMANN_J_PER_K
        * T_K
        / (AVOGADRO_PER_MOL * float(density_kg_per_m3) * mu**2)
    )


def fit_molten_pvdf_permittivity(
    csv_path: str | Path,
    *,
    epsilon_inf: float = 2.2,
    density_kg_per_m3: float = 1680.0,
    dipole_moment_D: float = 0.923,
    molar_mass_kg_per_mol: float = 0.064,
) -> KirkwoodFrohlichMeltFit:
    """Fit the Rui-2021 Kirkwood-Frohlich g(T) line from molten-PVDF epsilon(T).

    The function name is retained for API compatibility with v0.1.7.  The fitted
    quantity is now correctly g(T), not epsilon(T).
    """
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
        raise ValueError("need at least three molten-PVDF source points")

    T_C = np.asarray([p.temperature_C for p in selected], dtype=float)
    eps = np.asarray([p.value for p in selected], dtype=float)
    if np.unique(T_C).size != T_C.size:
        raise ValueError("duplicate molten-PVDF temperatures")
    if np.any(eps <= 0.0):
        raise ValueError("permittivity must be positive")

    g = np.asarray([
        _g_from_epsilon(
            e,
            T,
            epsilon_inf=epsilon_inf,
            density_kg_per_m3=density_kg_per_m3,
            dipole_moment_D=dipole_moment_D,
            molar_mass_kg_per_mol=molar_mass_kg_per_mol,
        )
        for T, e in zip(T_C, eps)
    ])
    T_K = T_C + 273.15
    slope, intercept = np.polyfit(T_K, g, 1)
    pred = slope * T_K + intercept
    residual = g - pred
    ss_res = float(np.sum(residual**2))
    ss_tot = float(np.sum((g - np.mean(g))**2))
    r_squared = 1.0 if ss_tot == 0.0 else 1.0 - ss_res / ss_tot
    rmse = float(np.sqrt(np.mean(residual**2)))

    return KirkwoodFrohlichMeltFit(
        epsilon_inf=float(epsilon_inf),
        density_kg_per_m3=float(density_kg_per_m3),
        dipole_moment_D=float(dipole_moment_D),
        molar_mass_kg_per_mol=float(molar_mass_kg_per_mol),
        g_slope_per_K=float(slope),
        g_intercept=float(intercept),
        g_r_squared=float(r_squared),
        g_rmse=rmse,
        temperature_min_C=float(T_C.min()),
        temperature_max_C=float(T_C.max()),
        n_points=int(T_C.size),
    )


def epsilon_iaf_from_melt_extrapolation(
    csv_path: str | Path,
    *,
    temperature_C: float,
    allowed_low_C: float = -30.0,
    allowed_high_C: float = 40.0,
) -> tuple[float, KirkwoodFrohlichMeltFit]:
    """Return the source-model IAF/amorphous epsilon(T) in the project BDS window."""
    T = float(temperature_C)
    if not (allowed_low_C <= T <= allowed_high_C):
        raise ValueError(
            f"epsilon_IAF source-model extrapolation is restricted to "
            f"[{allowed_low_C:g}, {allowed_high_C:g}] C"
        )
    fit = fit_molten_pvdf_permittivity(csv_path)
    return fit.epsilon(T), fit
