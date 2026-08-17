"""Source-faithful first reproduction step for Guo et al. (Nat. Commun. 2024).

Primary source
--------------
Mengfan Guo et al., "Electrically and mechanically driven rotation of polar
spirals in a relaxor ferroelectric polymer", Nature Communications 15, 348
(2024), DOI 10.1038/s41467-023-44395-5.

Supplementary Tables S2 and S3 report strong- and weak-anisotropy Landau
coefficients. The complete paper uses a vector sixth-order Landau polynomial,
but the expanded cross-term convention is not printed in a form we can safely
reconstruct from the table alone. This module therefore starts with a strict
one-axis slice

    f(P) = alpha1 P^2 + alpha11 P^4 + alpha111 P^6,

for which every cross term vanishes. No unreported coefficient is introduced.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class AxisLandauParameters:
    """Directly reported coefficients sufficient for a one-axis Landau slice."""

    label: str
    alpha1_prefactor_J_m_C2_K: float
    curie_temperature_C: float
    alpha11_J_m5_C4: float
    alpha111_J_m9_C6: float
    source_table: str

    def alpha1(self, temperature_C: float) -> float:
        """Return alpha1(T) in J m C^-2 using the source-reported law."""
        return float(self.alpha1_prefactor_J_m_C2_K) * (
            float(temperature_C) - float(self.curie_temperature_C)
        )


# Direct transcription from Guo et al. 2024 Supplementary Tables S2/S3.
# The alpha1 entry is printed as 1.412 (T - 42) x 10^5 J m C^-2.
GUO2024_STRONG_ANISOTROPY_AXIS = AxisLandauParameters(
    label="strong_anisotropy_S2",
    alpha1_prefactor_J_m_C2_K=1.412e5,
    curie_temperature_C=42.0,
    alpha11_J_m5_C4=-1.842e9,
    alpha111_J_m9_C6=2.585e11,
    source_table="Supplementary Table S2",
)

GUO2024_WEAK_ANISOTROPY_AXIS = AxisLandauParameters(
    label="weak_anisotropy_S3",
    alpha1_prefactor_J_m_C2_K=1.412e5,
    curie_temperature_C=42.0,
    alpha11_J_m5_C4=-1.842e8,
    alpha111_J_m9_C6=2.585e12,
    source_table="Supplementary Table S3",
)


def landau_axis_density(
    polarization_C_m2: np.ndarray | float,
    temperature_C: float,
    parameters: AxisLandauParameters,
) -> np.ndarray | float:
    """Evaluate the source-supported one-axis Landau density in J m^-3."""
    p = np.asarray(polarization_C_m2, dtype=float)
    a1 = parameters.alpha1(temperature_C)
    f = (
        a1 * p**2
        + parameters.alpha11_J_m5_C4 * p**4
        + parameters.alpha111_J_m9_C6 * p**6
    )
    if np.ndim(polarization_C_m2) == 0:
        return float(f)
    return f


def landau_axis_derivative(
    polarization_C_m2: np.ndarray | float,
    temperature_C: float,
    parameters: AxisLandauParameters,
) -> np.ndarray | float:
    """Return df/dP for the one-axis source-supported polynomial."""
    p = np.asarray(polarization_C_m2, dtype=float)
    a1 = parameters.alpha1(temperature_C)
    df = (
        2.0 * a1 * p
        + 4.0 * parameters.alpha11_J_m5_C4 * p**3
        + 6.0 * parameters.alpha111_J_m9_C6 * p**5
    )
    if np.ndim(polarization_C_m2) == 0:
        return float(df)
    return df


def axis_equilibrium(
    temperature_C: float,
    parameters: AxisLandauParameters,
) -> dict[str, float]:
    """Derive the stable non-zero one-axis minimum when it exists.

    This is a repository-derived consequence of the source coefficients, not a
    value quoted by Guo et al. For P != 0, stationarity gives

        3 alpha111 y^2 + 2 alpha11 y + alpha1 = 0,  y = P^2.

    All positive stationary roots are checked and only positive-curvature minima
    compete with P=0 for the global one-axis minimum.
    """
    a1 = parameters.alpha1(temperature_C)
    a11 = float(parameters.alpha11_J_m5_C4)
    a111 = float(parameters.alpha111_J_m9_C6)
    if not all(math.isfinite(v) for v in (a1, a11, a111)) or a111 <= 0.0:
        raise ValueError("Landau coefficients must be finite and alpha111 positive")

    discriminant = (2.0 * a11) ** 2 - 12.0 * a111 * a1
    candidates: list[tuple[float, float, float]] = [(0.0, 0.0, 2.0 * a1)]
    if discriminant >= 0.0:
        root_disc = math.sqrt(discriminant)
        for y in (
            (-2.0 * a11 + root_disc) / (6.0 * a111),
            (-2.0 * a11 - root_disc) / (6.0 * a111),
        ):
            if y <= 0.0:
                continue
            p = math.sqrt(y)
            f = float(landau_axis_density(p, temperature_C, parameters))
            curvature = float(
                2.0 * a1 + 12.0 * a11 * y + 30.0 * a111 * y**2
            )
            if curvature > 0.0:
                candidates.append((p, f, curvature))

    p_abs, f_min, curvature = min(candidates, key=lambda item: item[1])
    return {
        "temperature_C": float(temperature_C),
        "alpha1_J_m_C2": float(a1),
        "P_abs_min_C_m2": float(p_abs),
        "f_min_J_m3": float(f_min),
        # d^2 f / dP^2 has the same units as alpha1: J m C^-2.
        "curvature_at_min_J_m_C2": float(curvature),
    }
