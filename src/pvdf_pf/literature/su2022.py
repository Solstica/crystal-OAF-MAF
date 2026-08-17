"""Source-faithful homogeneous beta-PVDF checkpoint from Su et al. (2022).

Primary source:
    Yuanjie Su et al., Nature Communications 13, 4867 (2022),
    DOI 10.1038/s41467-022-32518-3.

The paper explicitly prints the PVDF bulk free-energy density

    f_bulk = alpha1 Px^2 + alpha2 Py^2 + alpha3 Pz^2
             + alpha33 Pz^4 + alpha333 Pz^6,

with a uniaxial ferroelectric z direction and paraelectric x/y directions.
Supplementary Table 3 provides all coefficients used below. No gradient or
kinetic coefficient is inferred here because those numerical values are not
reported in the provided source files.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class Su2022PVDFParameters:
    alpha1_J_m_C2: float
    alpha2_J_m_C2: float
    alpha3_prefactor_J_m_C2_K: float
    alpha3_zero_K: float
    alpha33_J_m5_C4: float
    alpha333_J_m9_C6: float
    Q11_m4_C2: float
    Q12_m4_C2: float
    Q44_m4_C2: float
    s11_m2_N: float
    s12_m2_N: float
    s44_m2_N: float
    source_table: str = "Supplementary Table 3"

    def alpha3(self, temperature_K: float) -> float:
        """Temperature-dependent z-axis quadratic coefficient."""
        temperature_K = float(temperature_K)
        if not math.isfinite(temperature_K):
            raise ValueError("temperature_K must be finite")
        return self.alpha3_prefactor_J_m_C2_K * (
            temperature_K - self.alpha3_zero_K
        )


SU2022_PVDF = Su2022PVDFParameters(
    alpha1_J_m_C2=5.647e9,
    alpha2_J_m_C2=5.647e9,
    alpha3_prefactor_J_m_C2_K=1.412e7,
    alpha3_zero_K=315.0,
    alpha33_J_m5_C4=-1.842e11,
    alpha333_J_m9_C6=2.585e13,
    Q11_m4_C2=-8.5,
    Q12_m4_C2=0.0,
    Q44_m4_C2=0.0,
    s11_m2_N=4.0e-10,
    s12_m2_N=1.11e-9,
    s44_m2_N=1.25e-9,
)


def pvdf_bulk_density(
    polarization_C_m2: np.ndarray,
    temperature_K: float,
    parameters: Su2022PVDFParameters = SU2022_PVDF,
) -> np.ndarray | float:
    """Evaluate Su et al. Eq. (4) for PVDF; returned units are J m^-3."""
    p = np.asarray(polarization_C_m2, dtype=float)
    if p.shape[-1] != 3:
        raise ValueError("polarization_C_m2 must have final dimension length 3")
    px, py, pz = np.moveaxis(p, -1, 0)
    f = (
        parameters.alpha1_J_m_C2 * px**2
        + parameters.alpha2_J_m_C2 * py**2
        + parameters.alpha3(temperature_K) * pz**2
        + parameters.alpha33_J_m5_C4 * pz**4
        + parameters.alpha333_J_m9_C6 * pz**6
    )
    if p.ndim == 1:
        return float(f)
    return f


def z_axis_stationary_points(
    temperature_K: float,
    parameters: Su2022PVDFParameters = SU2022_PVDF,
) -> list[dict[str, float | str]]:
    """Return nonnegative-Pz stationary points of the printed homogeneous model.

    Values are project-derived algebraic checkpoints at a declared temperature,
    not quantities quoted by Su et al.
    """
    a2 = parameters.alpha3(temperature_K)
    a4 = parameters.alpha33_J_m5_C4
    a6 = parameters.alpha333_J_m9_C6
    points: list[dict[str, float | str]] = []

    def add_point(p: float) -> None:
        p2 = p * p
        energy = a2 * p2 + a4 * p2**2 + a6 * p2**3
        curvature = 2.0 * a2 + 12.0 * a4 * p2 + 30.0 * a6 * p2**2
        if curvature > 0.0:
            kind = "minimum"
        elif curvature < 0.0:
            kind = "maximum"
        else:
            kind = "flat"
        points.append(
            {
                "Pz_C_m2": float(p),
                "f_J_m3": float(energy),
                "curvature_J_m_C2": float(curvature),
                "kind": kind,
            }
        )

    add_point(0.0)
    disc = (2.0 * a4) ** 2 - 12.0 * a6 * a2
    if disc >= 0.0:
        root_disc = math.sqrt(disc)
        for p2 in (
            (-2.0 * a4 - root_disc) / (6.0 * a6),
            (-2.0 * a4 + root_disc) / (6.0 * a6),
        ):
            if p2 > 0.0:
                add_point(math.sqrt(p2))
    return sorted(points, key=lambda item: float(item["Pz_C_m2"]))


def z_axis_global_minimum(
    temperature_K: float,
    parameters: Su2022PVDFParameters = SU2022_PVDF,
) -> dict[str, float | str]:
    """Return the lowest-energy nonnegative stationary state."""
    points = z_axis_stationary_points(temperature_K, parameters)
    return min(points, key=lambda item: float(item["f_J_m3"]))


def electrostatic_displacement(
    electric_field_V_m: np.ndarray,
    polarization_C_m2: np.ndarray,
    eps_b: float,
    eps0_F_m: float = 8.8541878128e-12,
) -> np.ndarray:
    """Evaluate D = eps0*eps_b*E + P from Su et al. Eq. (5).

    The source states the equilibrium equation div(D)=0. This helper deliberately
    requires eps_b as an explicit argument because the paper varies it with MXene
    loading; the repository does not invent a default value.
    """
    e = np.asarray(electric_field_V_m, dtype=float)
    p = np.asarray(polarization_C_m2, dtype=float)
    if e.shape != p.shape or e.shape[-1] != 3:
        raise ValueError("E and P must have identical shapes ending in length 3")
    if not math.isfinite(float(eps_b)) or float(eps_b) <= 0.0:
        raise ValueError("eps_b must be finite and positive")
    return float(eps0_F_m) * float(eps_b) * e + p
