"""Dimensionless-parameter audit for Ahluwalia et al. PRB 78, 054110 (2008).

This module reconstructs the nondimensional quantities printed around Eqs. (11)-(12)
from the paper's Table-II LGD coefficients and the reported 300 K noise match.

Important provenance distinction:
- the algebraic rescaling relations are SOURCE_EQUATIONS;
- the reported noise amplitudes and physical grid/cell sizes are DIRECT_SOURCE_VALUES;
- the characteristic length inferred by inverting the reported z-noise relation is
  PROJECT_INFERRED_FROM_SOURCE_EQ11_AND_MATCHED_NOISE.

The inferred length is not silently identified with the 2.16 nm continuum grid spacing.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import math

from pvdf_pf.literature.ahluwalia2008 import (
    GRADIENT_CHECKPOINT,
    KINETIC_CHECKPOINT,
    TABLE_II_PRINTED,
    ferroelectric_stationary_magnitude,
    kinetic_ratios,
    quadratic_coefficient,
)

K_B_J_K = 1.380649e-23
EPS0_F_M = 8.8541878128e-12


@dataclass(frozen=True)
class AhluwaliaDimensionlessCheckpoint:
    temperature_K: float
    polarization_scale_C_m2: float
    equilibrium_Pr_C_m2: float
    equilibrium_w_dimensionless: float
    initial_w_dimensionless: float
    alpha_prime: float
    transverse_susceptibility_m_F: float
    alpha_xx_prime: float
    alpha_yy_prime: float
    m_Gamma_x_over_Gamma_z: float
    n_Gamma_y_over_Gamma_z: float
    E_prime: float
    inferred_length_scale_m: float
    inferred_length_scale_nm: float
    K1_prime: float
    K2_prime: float
    K3_prime: float
    physical_grid_spacing_nm: float
    grid_spacing_dimensionless: float
    physical_cell_length_nm: float
    cell_length_dimensionless: float
    noise_x_dimensionless: float
    noise_y_dimensionless: float
    noise_z_dimensionless: float
    ps_per_tstar: float


def polarization_scale(parameters=TABLE_II_PRINTED) -> float:
    """Return sqrt(beta/gamma), the source polarization scaling relation."""
    return float(math.sqrt(parameters.beta_J_m5_C4 / parameters.gamma_J_m9_C6))


def inferred_length_from_noise_z(
    temperature_K: float,
    noise_z_dimensionless: float = KINETIC_CHECKPOINT.noise_z_dimensionless,
    parameters=TABLE_II_PRINTED,
) -> float:
    """Invert the source z-noise rescaling to obtain the characteristic length.

    Around Eq. (11), the source gives
        eps_z_tilde = sqrt(2 k_B T / (beta P_scale^4 xi^3)).
    Hence
        xi = [2 k_B T / (beta P_scale^4 eps_z_tilde^2)]^(1/3).

    The returned value is a project-derived inference from a printed equation and a
    reported matched noise amplitude, not a separately quoted length in the paper.
    """
    t = float(temperature_K)
    eta = float(noise_z_dimensionless)
    if not math.isfinite(t) or t <= 0.0:
        raise ValueError("temperature_K must be finite and positive")
    if not math.isfinite(eta) or eta <= 0.0:
        raise ValueError("noise_z_dimensionless must be finite and positive")
    pscale = polarization_scale(parameters)
    return float(
        (
            2.0
            * K_B_J_K
            * t
            / (parameters.beta_J_m5_C4 * pscale**4 * eta**2)
        )
        ** (1.0 / 3.0)
    )


def dimensionless_checkpoint_300K() -> AhluwaliaDimensionlessCheckpoint:
    """Reconstruct the source's 300 K nondimensional coefficient set.

    The paper assumes chi_xx=chi_yy=chi_zz, with chi_zz evaluated as the inverse
    homogeneous Landau curvature at the remnant polarization. The present audit uses
    the stable homogeneous 300 K branch from the directly printed Table-II coefficients.
    """
    t = 300.0
    params = TABLE_II_PRINTED
    pscale = polarization_scale(params)
    equilibrium = ferroelectric_stationary_magnitude(t, params)
    if equilibrium is None:
        raise ValueError("expected a ferroelectric stationary state at 300 K")
    pr = float(equilibrium["P_C_m2"])
    curvature = float(equilibrium["curvature_J_m_C2"])
    chi = 1.0 / curvature

    beta_ps2 = params.beta_J_m5_C4 * pscale**2
    alpha_prime = quadratic_coefficient(t, params) / beta_ps2
    alpha_transverse_prime = 1.0 / (beta_ps2 * chi)
    e_prime = 1.0 / (EPS0_F_M * beta_ps2)

    xi = inferred_length_from_noise_z(t, KINETIC_CHECKPOINT.noise_z_dimensionless, params)
    kprime = GRADIENT_CHECKPOINT.K1_J_m3_C2 / (beta_ps2 * xi**2)
    ratios = kinetic_ratios()

    grid_nm = KINETIC_CHECKPOINT.tdgl_grid_length_nm
    cell_nm = KINETIC_CHECKPOINT.tdgl_cell_length_nm
    return AhluwaliaDimensionlessCheckpoint(
        temperature_K=t,
        polarization_scale_C_m2=pscale,
        equilibrium_Pr_C_m2=pr,
        equilibrium_w_dimensionless=pr / pscale,
        initial_w_dimensionless=KINETIC_CHECKPOINT.initial_Pz_C_m2 / pscale,
        alpha_prime=float(alpha_prime),
        transverse_susceptibility_m_F=float(chi),
        alpha_xx_prime=float(alpha_transverse_prime),
        alpha_yy_prime=float(alpha_transverse_prime),
        m_Gamma_x_over_Gamma_z=float(ratios["m_Gamma_x_over_Gamma_z"]),
        n_Gamma_y_over_Gamma_z=float(ratios["n_Gamma_y_over_Gamma_z"]),
        E_prime=float(e_prime),
        inferred_length_scale_m=float(xi),
        inferred_length_scale_nm=float(xi * 1e9),
        K1_prime=float(kprime),
        K2_prime=float(kprime),
        K3_prime=float(kprime),
        physical_grid_spacing_nm=float(grid_nm),
        grid_spacing_dimensionless=float(grid_nm * 1e-9 / xi),
        physical_cell_length_nm=float(cell_nm),
        cell_length_dimensionless=float(cell_nm * 1e-9 / xi),
        noise_x_dimensionless=KINETIC_CHECKPOINT.noise_x_dimensionless,
        noise_y_dimensionless=KINETIC_CHECKPOINT.noise_y_dimensionless,
        noise_z_dimensionless=KINETIC_CHECKPOINT.noise_z_dimensionless,
        ps_per_tstar=float(ratios["ps_per_tstar_from_equilibration_match"]),
    )


def dimensionless_checkpoint_report() -> dict[str, object]:
    checkpoint = dimensionless_checkpoint_300K()
    return {
        "checkpoint": asdict(checkpoint),
        "provenance": {
            "polarization_scale": "SOURCE_EQ11_RESCALING_RELATION_USING_PRINTED_TABLE_II",
            "noise_amplitudes": "DIRECT_SOURCE_300K_MATCH",
            "physical_grid_and_cell": "DIRECT_SOURCE",
            "inferred_length_scale": "PROJECT_INFERRED_FROM_SOURCE_EQ11_AND_MATCHED_NOISE",
            "dimensionless_coefficients": "PROJECT_DERIVED_FROM_SOURCE_EQUATIONS_AND_DIRECT_VALUES",
            "K3": GRADIENT_CHECKPOINT.K3_status,
        },
    }
