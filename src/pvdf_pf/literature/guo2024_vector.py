"""Source-constrained reconstruction of the Guo et al. (2024) homogeneous
three-component Landau polynomial used for Supplementary Fig. S16.

Guo et al. print the tensor form of the sixth-order Landau energy and tabulate
the contracted coefficients in Supplementary Tables S2/S3, but do not print the
expanded polynomial. The expansion used here is the conventional cubic
sixth-order Landau-Devonshire form

    a1 * sum(P_i^2)
  + a11 * sum(P_i^4)
  + a12 * sum_{i<j}(P_i^2 P_j^2)
  + a111 * sum(P_i^6)
  + a112 * sum_i P_i^4 * sum_{j != i} P_j^2
  + a123 * P_x^2 P_y^2 P_z^2.

This coefficient convention is cross-checked against the explicitly expanded
sixth-order phase-field polynomial printed by Su et al., Nature Communications
13, 4867 (2022), Eq. (3). It is therefore tagged as a published-convention
cross-check, not as a verbatim expanded equation from Guo et al.

No fitting or project-defined anisotropy coefficient is introduced.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class CubicSixthOrderLandauParameters:
    label: str
    alpha1_prefactor_J_m_C2_K: float
    curie_temperature_C: float
    alpha11_J_m5_C4: float
    alpha12_J_m5_C4: float
    alpha111_J_m9_C6: float
    alpha112_J_m9_C6: float
    alpha123_J_m9_C6: float
    source_table: str

    def alpha1(self, temperature_C: float) -> float:
        return float(self.alpha1_prefactor_J_m_C2_K) * (
            float(temperature_C) - float(self.curie_temperature_C)
        )


GUO2024_STRONG_ANISOTROPY_VECTOR = CubicSixthOrderLandauParameters(
    label="strong_anisotropy_S2",
    alpha1_prefactor_J_m_C2_K=1.412e5,
    curie_temperature_C=42.0,
    alpha11_J_m5_C4=-1.842e9,
    alpha12_J_m5_C4=0.0,
    alpha111_J_m9_C6=2.585e11,
    alpha112_J_m9_C6=0.0,
    alpha123_J_m9_C6=0.0,
    source_table="Supplementary Table S2",
)

GUO2024_WEAK_ANISOTROPY_VECTOR = CubicSixthOrderLandauParameters(
    label="weak_anisotropy_S3",
    alpha1_prefactor_J_m_C2_K=1.412e5,
    curie_temperature_C=42.0,
    alpha11_J_m5_C4=-1.842e8,
    alpha12_J_m5_C4=-1.4736e9,
    alpha111_J_m9_C6=2.585e12,
    alpha112_J_m9_C6=9.6e12,
    alpha123_J_m9_C6=1.0857e13,
    source_table="Supplementary Table S3",
)


def vector_landau_density(
    polarization_C_m2: np.ndarray,
    temperature_C: float,
    parameters: CubicSixthOrderLandauParameters,
) -> np.ndarray | float:
    """Evaluate the conventional expanded sixth-order Landau density.

    The final dimension of ``polarization_C_m2`` must have length 3 and store
    (Px, Py, Pz). Returned units are J m^-3.
    """
    p = np.asarray(polarization_C_m2, dtype=float)
    if p.shape[-1] != 3:
        raise ValueError("polarization_C_m2 must have final dimension length 3")
    x, y, z = np.moveaxis(p, -1, 0)
    x2, y2, z2 = x * x, y * y, z * z
    a1 = parameters.alpha1(temperature_C)
    f = (
        a1 * (x2 + y2 + z2)
        + parameters.alpha11_J_m5_C4 * (x2**2 + y2**2 + z2**2)
        + parameters.alpha12_J_m5_C4 * (x2 * y2 + x2 * z2 + y2 * z2)
        + parameters.alpha111_J_m9_C6 * (x2**3 + y2**3 + z2**3)
        + parameters.alpha112_J_m9_C6
        * (
            x2**2 * (y2 + z2)
            + y2**2 * (z2 + x2)
            + z2**2 * (x2 + y2)
        )
        + parameters.alpha123_J_m9_C6 * x2 * y2 * z2
    )
    if p.ndim == 1:
        return float(f)
    return f


def directional_polynomial_coefficients(
    direction: np.ndarray,
    temperature_C: float,
    parameters: CubicSixthOrderLandauParameters,
) -> tuple[float, float, float]:
    """Return A2,A4,A6 for f(r*n)=A2 r^2 + A4 r^4 + A6 r^6."""
    n = np.asarray(direction, dtype=float)
    if n.shape != (3,):
        raise ValueError("direction must have shape (3,)")
    norm = float(np.linalg.norm(n))
    if not math.isfinite(norm) or norm <= 0.0:
        raise ValueError("direction must be finite and non-zero")
    x, y, z = n / norm
    x2, y2, z2 = x * x, y * y, z * z
    s4 = x2**2 + y2**2 + z2**2
    pair4 = x2 * y2 + x2 * z2 + y2 * z2
    s6 = x2**3 + y2**3 + z2**3
    mixed6 = x2**2 * (y2 + z2) + y2**2 * (z2 + x2) + z2**2 * (x2 + y2)
    triple6 = x2 * y2 * z2
    return (
        float(parameters.alpha1(temperature_C)),
        float(parameters.alpha11_J_m5_C4 * s4 + parameters.alpha12_J_m5_C4 * pair4),
        float(
            parameters.alpha111_J_m9_C6 * s6
            + parameters.alpha112_J_m9_C6 * mixed6
            + parameters.alpha123_J_m9_C6 * triple6
        ),
    )


def equilibrium_along_direction(
    direction: np.ndarray,
    temperature_C: float,
    parameters: CubicSixthOrderLandauParameters,
) -> dict[str, float]:
    """Global homogeneous minimum constrained to a fixed polarization direction.

    This is a deterministic project-derived checkpoint from source coefficients.
    It is not a value quoted by Guo et al.
    """
    a2, a4, a6 = directional_polynomial_coefficients(
        direction, temperature_C, parameters
    )
    if not all(math.isfinite(v) for v in (a2, a4, a6)) or a6 <= 0.0:
        raise ValueError("directional sixth-order coefficient must be positive")
    disc = (2.0 * a4) ** 2 - 12.0 * a6 * a2
    candidates: list[tuple[float, float, float]] = [(0.0, 0.0, 2.0 * a2)]
    if disc >= 0.0:
        root_disc = math.sqrt(disc)
        for p2 in (
            (-2.0 * a4 + root_disc) / (6.0 * a6),
            (-2.0 * a4 - root_disc) / (6.0 * a6),
        ):
            if p2 <= 0.0:
                continue
            r = math.sqrt(p2)
            energy = a2 * p2 + a4 * p2**2 + a6 * p2**3
            curvature = 2.0 * a2 + 12.0 * a4 * p2 + 30.0 * a6 * p2**2
            if curvature > 0.0:
                candidates.append((r, energy, curvature))
    r, energy, curvature = min(candidates, key=lambda item: item[1])
    return {
        "temperature_C": float(temperature_C),
        "P_abs_min_C_m2": float(r),
        "f_min_J_m3": float(energy),
        "radial_curvature_J_m_C2": float(curvature),
        "A2_J_m_C2": float(a2),
        "A4_J_m5_C4": float(a4),
        "A6_J_m9_C6": float(a6),
    }


def fibonacci_sphere(n_points: int) -> np.ndarray:
    """Deterministic approximately uniform unit directions."""
    if n_points < 2:
        raise ValueError("n_points must be >= 2")
    i = np.arange(n_points, dtype=float)
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    z = 1.0 - 2.0 * (i + 0.5) / float(n_points)
    radius_xy = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    phi = golden_angle * i
    return np.column_stack(
        (radius_xy * np.cos(phi), radius_xy * np.sin(phi), z)
    )


def directional_anisotropy_summary(
    temperature_C: float,
    parameters: CubicSixthOrderLandauParameters,
    n_points: int = 4096,
) -> dict[str, float]:
    """Summarize radial-minimum and minimum-depth angular anisotropy."""
    directions = fibonacci_sphere(n_points)
    radii = np.empty(n_points, dtype=float)
    energies = np.empty(n_points, dtype=float)
    for idx, direction in enumerate(directions):
        result = equilibrium_along_direction(direction, temperature_C, parameters)
        radii[idx] = result["P_abs_min_C_m2"]
        energies[idx] = result["f_min_J_m3"]

    depths = -energies
    if np.any(radii <= 0.0) or np.any(depths <= 0.0):
        raise ValueError("expected ferroelectric directional minima at this temperature")
    return {
        "n_points": int(n_points),
        "P_radius_min_C_m2": float(radii.min()),
        "P_radius_max_C_m2": float(radii.max()),
        "P_radius_max_over_min": float(radii.max() / radii.min()),
        "P_radius_range_over_mean": float((radii.max() - radii.min()) / radii.mean()),
        "well_depth_min_J_m3": float(depths.min()),
        "well_depth_max_J_m3": float(depths.max()),
        "well_depth_max_over_min": float(depths.max() / depths.min()),
        "well_depth_range_over_mean": float(
            (depths.max() - depths.min()) / depths.mean()
        ),
    }
