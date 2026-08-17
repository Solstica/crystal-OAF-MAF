"""Source-faithful checkpoints for Ahluwalia et al., Phys. Rev. B 78, 054110 (2008).

The paper is used here as the multiscale-method anchor because molecular-dynamics
outputs are explicitly transferred into an LGD/TDGL continuum model. Directly
reported values are kept separate from values re-derived from the rounded printed
tables so that source precision is not silently upgraded.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math


@dataclass(frozen=True)
class AhluwaliaTableI:
    P0_C_m2: float
    Tc_K: float
    Pc_C_m2: float
    f0_J_m3: float


@dataclass(frozen=True)
class AhluwaliaLGDParameters:
    alpha0_J_m_C2_K: float
    T0_K: float
    beta_J_m5_C4: float
    gamma_J_m9_C6: float


TABLE_I_PRINTED = AhluwaliaTableI(
    P0_C_m2=0.111,
    Tc_K=450.0,
    Pc_C_m2=0.088,
    f0_J_m3=-2.12e8,
)

TABLE_II_PRINTED = AhluwaliaLGDParameters(
    alpha0_J_m_C2_K=9.02e7,
    T0_K=252.0,
    beta_J_m5_C4=9.20e12,
    gamma_J_m9_C6=8.88e14,
)


@dataclass(frozen=True)
class AhluwaliaGradientCheckpoint:
    xi1_m: float
    xi2_m: float
    K1_J_m3_C2: float
    K2_J_m3_C2: float
    K3_J_m3_C2: float
    K3_status: str


GRADIENT_CHECKPOINT = AhluwaliaGradientCheckpoint(
    xi1_m=0.4e-9,
    xi2_m=0.4e-9,
    K1_J_m3_C2=2.108e-8,
    K2_J_m3_C2=2.108e-8,
    K3_J_m3_C2=2.108e-8,
    K3_status="SOURCE_COMPUTATIONAL_CONVENIENCE_SET_EQUAL_TO_K1_K2",
)


@dataclass(frozen=True)
class AhluwaliaKineticCheckpoint:
    noise_x_dimensionless: float
    noise_y_dimensionless: float
    noise_z_dimensionless: float
    tdgl_equilibration_tstar: float
    md_equilibration_ps: float
    reported_smallest_tdgl_time_ps: float
    tdgl_cell_length_nm: float
    tdgl_grid_length_nm: float
    initial_Pz_C_m2: float


KINETIC_CHECKPOINT = AhluwaliaKineticCheckpoint(
    noise_x_dimensionless=0.0711,
    noise_y_dimensionless=1.0669,
    noise_z_dimensionless=0.3556,
    tdgl_equilibration_tstar=5.0,
    md_equilibration_ps=45.0,
    reported_smallest_tdgl_time_ps=9.0,
    tdgl_cell_length_nm=138.24,
    tdgl_grid_length_nm=2.16,
    initial_Pz_C_m2=0.114,
)


def derive_lgd_from_table_i(
    table_i: AhluwaliaTableI = TABLE_I_PRINTED,
) -> AhluwaliaLGDParameters:
    """Apply Eqs. (2)-(3) to the printed Table-I values.

    The paper's Table I is rounded. Consequently this calculation is expected to
    differ slightly from the independently printed Table II and must not replace it.
    """
    P0 = float(table_i.P0_C_m2)
    Pc = float(table_i.Pc_C_m2)
    Tc = float(table_i.Tc_K)
    f0 = float(table_i.f0_J_m3)
    if not all(math.isfinite(v) for v in (P0, Pc, Tc, f0)):
        raise ValueError("Table-I values must be finite")
    if P0 <= 0.0 or Pc <= 0.0 or Tc <= 0.0 or f0 >= 0.0 or P0 <= Pc:
        raise ValueError("Table-I values are outside the source model regime")

    denominator = Pc**2 * P0**4 - P0**6
    beta = 4.0 * f0 * Pc**2 / denominator
    gamma = 3.0 * f0 / denominator
    delta = beta * (Pc**2 - P0**2) - gamma * (Pc**4 - P0**4)
    alpha0 = delta / Tc
    T0 = -Tc * (beta * P0**2 - gamma * P0**4) / delta
    return AhluwaliaLGDParameters(
        alpha0_J_m_C2_K=float(alpha0),
        T0_K=float(T0),
        beta_J_m5_C4=float(beta),
        gamma_J_m9_C6=float(gamma),
    )


def relative_difference(derived: float, printed: float) -> float:
    return float((derived - printed) / printed)


def table_ii_rounding_audit() -> dict[str, dict[str, float]]:
    derived = derive_lgd_from_table_i()
    result: dict[str, dict[str, float]] = {}
    for key, printed_value in asdict(TABLE_II_PRINTED).items():
        derived_value = getattr(derived, key)
        result[key] = {
            "printed": float(printed_value),
            "derived_from_rounded_table_I": float(derived_value),
            "relative_difference": relative_difference(derived_value, printed_value),
        }
    return result


def kinetic_ratios() -> dict[str, float]:
    """Derive kinetic-coefficient ratios encoded by the source noise amplitudes.

    Eq. (11) states noise_x = sqrt(m) noise_z and noise_y = sqrt(n) noise_z,
    where Gamma_x=m Gamma_z and Gamma_y=n Gamma_z.
    """
    k = KINETIC_CHECKPOINT
    m = (k.noise_x_dimensionless / k.noise_z_dimensionless) ** 2
    n = (k.noise_y_dimensionless / k.noise_z_dimensionless) ** 2
    tstar_to_ps = k.md_equilibration_ps / k.tdgl_equilibration_tstar
    return {
        "m_Gamma_x_over_Gamma_z": float(m),
        "n_Gamma_y_over_Gamma_z": float(n),
        "ps_per_tstar_from_equilibration_match": float(tstar_to_ps),
    }


def quadratic_coefficient(
    temperature_K: float,
    parameters: AhluwaliaLGDParameters = TABLE_II_PRINTED,
) -> float:
    """Return alpha0*(T-T0) in the source Eq. (1)."""
    temperature_K = float(temperature_K)
    if not math.isfinite(temperature_K):
        raise ValueError("temperature_K must be finite")
    return parameters.alpha0_J_m_C2_K * (temperature_K - parameters.T0_K)


def homogeneous_landau_density(
    polarization_C_m2: float,
    temperature_K: float,
    parameters: AhluwaliaLGDParameters = TABLE_II_PRINTED,
) -> float:
    """Evaluate source Eq. (1), f(P,T), in J m^-3."""
    p = float(polarization_C_m2)
    if not math.isfinite(p):
        raise ValueError("polarization_C_m2 must be finite")
    a = quadratic_coefficient(temperature_K, parameters)
    return float(
        0.5 * a * p**2
        - 0.25 * parameters.beta_J_m5_C4 * p**4
        + (1.0 / 6.0) * parameters.gamma_J_m9_C6 * p**6
    )


def homogeneous_equation_of_state(
    polarization_C_m2: float,
    temperature_K: float,
    parameters: AhluwaliaLGDParameters = TABLE_II_PRINTED,
) -> float:
    """Return E=d f/dP for homogeneous switching, in V m^-1."""
    p = float(polarization_C_m2)
    a = quadratic_coefficient(temperature_K, parameters)
    return float(
        a * p
        - parameters.beta_J_m5_C4 * p**3
        + parameters.gamma_J_m9_C6 * p**5
    )


def ferroelectric_stationary_magnitude(
    temperature_K: float,
    parameters: AhluwaliaLGDParameters = TABLE_II_PRINTED,
) -> dict[str, float] | None:
    """Positive nonzero stationary branch used for the source P(T) curve.

    This implements the positive-root analytical branch following Eq. (1). Above
    its spinodal the nonzero branch does not exist and ``None`` is returned.
    """
    a = quadratic_coefficient(temperature_K, parameters)
    beta = parameters.beta_J_m5_C4
    gamma = parameters.gamma_J_m9_C6
    disc = beta**2 - 4.0 * gamma * a
    if disc < 0.0:
        return None
    p2 = (beta + math.sqrt(disc)) / (2.0 * gamma)
    if p2 <= 0.0:
        return None
    p = math.sqrt(p2)
    curvature = a - 3.0 * beta * p2 + 5.0 * gamma * p2**2
    return {
        "temperature_K": float(temperature_K),
        "P_C_m2": float(p),
        "f_J_m3": homogeneous_landau_density(p, temperature_K, parameters),
        "curvature_J_m_C2": float(curvature),
    }


def transition_implied_by_printed_table_ii(
    parameters: AhluwaliaLGDParameters = TABLE_II_PRINTED,
) -> dict[str, float]:
    """First-order coexistence point implied by the rounded printed Table II.

    For f=a P^2/2-beta P^4/4+gamma P^6/6, coexistence with P=0 occurs at
    P_c^2=3 beta/(4 gamma) and a_c=3 beta^2/(16 gamma).
    """
    beta = parameters.beta_J_m5_C4
    gamma = parameters.gamma_J_m9_C6
    p_c = math.sqrt(3.0 * beta / (4.0 * gamma))
    a_c = 3.0 * beta**2 / (16.0 * gamma)
    t_c = parameters.T0_K + a_c / parameters.alpha0_J_m_C2_K
    return {
        "Tc_K": float(t_c),
        "Pc_C_m2": float(p_c),
        "f_at_Pc_J_m3": homogeneous_landau_density(p_c, t_c, parameters),
    }


def intrinsic_coercive_spinodal(
    temperature_K: float,
    parameters: AhluwaliaLGDParameters = TABLE_II_PRINTED,
) -> dict[str, float]:
    """Homogeneous positive-branch spinodal underlying the intrinsic P-E loop.

    The field-controlled equation of state is E(P)=df/dP. Loss of stability of
    the positive ferroelectric branch satisfies dE/dP=0. The returned signed
    switching field is negative for switching from +P toward -P.
    """
    a = quadratic_coefficient(temperature_K, parameters)
    beta = parameters.beta_J_m5_C4
    gamma = parameters.gamma_J_m9_C6
    disc = (3.0 * beta) ** 2 - 20.0 * gamma * a
    if disc < 0.0:
        raise ValueError("no real homogeneous spinodal at this temperature")
    roots = [
        (3.0 * beta - math.sqrt(disc)) / (10.0 * gamma),
        (3.0 * beta + math.sqrt(disc)) / (10.0 * gamma),
    ]
    positive = [q for q in roots if q > 0.0]
    if not positive:
        raise ValueError("no positive ferroelectric spinodal at this temperature")
    p2 = max(positive)
    p = math.sqrt(p2)
    field = homogeneous_equation_of_state(p, temperature_K, parameters)
    return {
        "temperature_K": float(temperature_K),
        "P_spinodal_C_m2": float(p),
        "E_switch_from_positive_V_m": float(field),
        "Ec_magnitude_V_m": float(abs(field)),
    }
