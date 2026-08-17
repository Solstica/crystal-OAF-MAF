"""Frozen ferroelectric polarization source for the v0.1.16 orientation audit.

The generalized-Debye polarization remains a physical-time auxiliary response.  A
separate prescribed P_FE,z(r) can be added to Gauss' law without evolving the TDGL
order parameter.  This isolates the electrostatic consequence of a remanent
polarization hypothesis before any switching kinetics are introduced.

For the Huang et al. BOPVDF lamellae, the beta-chain / lamellar normal is primarily
along MD whereas the macroscopic poled polarization is along the film normal.  In
the ideal aligned limit P_FE is therefore tangent to the crystal/OAF interfaces.
A scalar Pz that varies only along MD has div(P_FE)=0; any field redistribution from
P_FE then requires tilt, waviness, a z-component of the morphology normal, or other
spatial nonuniformity along the polarization direction.
"""
from __future__ import annotations

import numpy as np

from pvdf_pf.morphology.three_phase import CRYSTAL, OAF
from pvdf_pf.physics.coupled_local_field import (
    CoupledRelaxationStep,
    solve_periodic_local_field_scalar,
)
from pvdf_pf.physics.dipolar import EPS0
from pvdf_pf.physics.generalized_debye import GeneralizedDebyeBank


def lamellar_phase_coordinate(
    grid,
    *,
    winding_z: int,
    winding_x: int,
    offset: float = 0.0,
) -> np.ndarray:
    """Return the periodic phase coordinate used by periodic_winding_three_phase."""
    if winding_z == 0 and winding_x == 0:
        raise ValueError("at least one winding number must be non-zero")
    iz = np.arange(grid.nz, dtype=float)[:, None] / grid.nz
    ix = np.arange(grid.nx, dtype=float)[None, :] / grid.nx
    return (winding_z * iz + winding_x * ix - float(offset)) % 1.0


def build_lamellar_frozen_polarization_z(
    phase_map: np.ndarray,
    grid,
    *,
    crystal_fraction: float,
    oaf_fraction: float,
    p_beta_C_m2: float,
    p_oaf_mean_C_m2: float,
    winding_z: int,
    winding_x: int,
    oaf_profile: str = "uniform",
    decay_length_fraction_of_half_oaf: float = 0.35,
    offset: float = 0.0,
) -> np.ndarray:
    """Build a prescribed z-polarization on beta crystal and OAF cells.

    `oaf_profile="interface_decay"` makes the OAF magnitude largest at the
    crystal/OAF interface and exponentially weaker toward the IAF side.  The
    discrete OAF mean is rescaled exactly to `p_oaf_mean_C_m2`, so uniform and
    decaying profiles can be compared at fixed total OAF dipole moment.

    The profile shape is a project spatial hypothesis.  Huang et al.'s MD supports
    stronger crystal-proximal orientational constraint qualitatively, but its
    enlarged lateral unit cell and 2 ns production window do not calibrate an
    absolute decay length or local polarization magnitude.
    """
    phase = np.asarray(phase_map)
    if phase.shape != grid.shape:
        raise ValueError("phase_map must match grid.shape")
    fc = float(crystal_fraction)
    fo = float(oaf_fraction)
    if fc <= 0.0 or fo <= 0.0 or fc + fo >= 1.0:
        raise ValueError("crystal_fraction and oaf_fraction must leave a positive IAF fraction")
    if not np.isfinite(p_beta_C_m2) or not np.isfinite(p_oaf_mean_C_m2):
        raise ValueError("polarization magnitudes must be finite")

    P = np.zeros(grid.shape, dtype=float)
    crystal_mask = phase == CRYSTAL
    oaf_mask = phase == OAF
    if not np.any(crystal_mask) or not np.any(oaf_mask):
        raise ValueError("phase_map must contain beta crystal and OAF cells")
    P[crystal_mask] = float(p_beta_C_m2)

    profile = str(oaf_profile).lower()
    if profile == "uniform":
        P[oaf_mask] = float(p_oaf_mean_C_m2)
        return P
    if profile != "interface_decay":
        raise ValueError("oaf_profile must be 'uniform' or 'interface_decay'")

    lam = float(decay_length_fraction_of_half_oaf)
    if not np.isfinite(lam) or lam <= 0.0:
        raise ValueError("decay_length_fraction_of_half_oaf must be finite and positive")

    u = lamellar_phase_coordinate(
        grid,
        winding_z=winding_z,
        winding_x=winding_x,
        offset=offset,
    )
    half_oaf = 0.5 * fo
    c0 = half_oaf
    c1 = c0 + fc

    distance = np.full(grid.shape, np.nan, dtype=float)
    left = oaf_mask & (u < c0)
    right = oaf_mask & (u >= c1)
    distance[left] = (c0 - u[left]) / half_oaf
    distance[right] = (u[right] - c1) / half_oaf
    if np.any(~np.isfinite(distance[oaf_mask])):
        raise ValueError("phase_map is inconsistent with the supplied laminate fractions/winding")

    raw = np.exp(-distance[oaf_mask] / lam)
    raw_mean = float(np.mean(raw))
    if raw_mean <= 0.0 or not np.isfinite(raw_mean):
        raise RuntimeError("invalid OAF decay profile normalization")
    P[oaf_mask] = float(p_oaf_mean_C_m2) * raw / raw_mean
    return P


def advance_generalized_debye_with_frozen_source(
    P_modes_old: np.ndarray,
    eps_background_r: np.ndarray,
    amorphous_mask: np.ndarray,
    P_frozen_z: np.ndarray,
    grid,
    *,
    E_external_z: float,
    dt_s: float,
    bank: GeneralizedDebyeBank,
    eps0: float = EPS0,
    tol: float = 1e-10,
    maxiter: int = 4000,
) -> CoupledRelaxationStep:
    """Advance one Debye step with a prescribed, non-evolving P_FE,z source.

    The algebraic ZOH elimination is identical to v0.1.15.  Only the explicit
    polarization entering Gauss' law changes from P_memory to
    P_memory + P_frozen_z.  P_frozen_z is then included when reconstructing D.
    """
    bank.validate()
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive")

    eps_b = np.asarray(eps_background_r, dtype=float)
    mask = np.asarray(amorphous_mask, dtype=bool)
    P_frozen = np.asarray(P_frozen_z, dtype=float)
    if eps_b.shape != grid.shape or mask.shape != grid.shape or P_frozen.shape != grid.shape:
        raise ValueError("eps_background_r, amorphous_mask and P_frozen_z must match grid.shape")
    if np.any(eps_b <= 0.0) or np.any(~np.isfinite(eps_b)):
        raise ValueError("eps_background_r must be finite and positive")
    if np.any(~np.isfinite(P_frozen)):
        raise ValueError("P_frozen_z must be finite")

    P_old = np.asarray(P_modes_old, dtype=float)
    expected = (bank.n_modes,) + grid.shape
    if P_old.shape != expected:
        raise ValueError(f"P_modes_old shape must be {expected}, got {P_old.shape}")
    if np.any(~np.isfinite(P_old)):
        raise ValueError("P_modes_old must be finite")

    tau = np.asarray(bank.tau_modes_s, dtype=float).reshape((-1, 1, 1))
    strength = np.asarray(bank.delta_epsilon_modes, dtype=float).reshape((-1, 1, 1))
    a = np.exp(-float(dt_s) / tau)
    b = 1.0 - a
    mask3 = mask[None, :, :]
    P_old_masked = np.where(mask3, P_old, 0.0)

    memory_modes = a * P_old_masked
    P_memory = np.sum(memory_modes, axis=0)
    step_increment = float(np.sum(b[:, 0, 0] * strength[:, 0, 0]))
    eps_step = eps_b + mask.astype(float) * step_increment

    field = solve_periodic_local_field_scalar(
        P_memory + P_frozen,
        eps_step,
        grid,
        E_external_z=float(E_external_z),
        eps0=eps0,
        tol=tol,
        maxiter=maxiter,
    )

    driven = b * float(eps0) * strength * field.E_z[None, :, :]
    P_new = np.where(mask3, memory_modes + driven, 0.0)
    P_total = np.sum(P_new, axis=0)

    D_cell = float(eps0) * eps_b * field.E_z + P_total + P_frozen
    D_mean_cell = float(np.mean(D_cell))
    scale = max(abs(field.D_z_mean_face), abs(D_mean_cell), 1e-30)
    mismatch = float(abs(field.D_z_mean_face - D_mean_cell) / scale)

    return CoupledRelaxationStep(
        P_modes=P_new,
        P_total=P_total,
        potential=field.potential,
        E_z=field.E_z,
        E_x=field.E_x,
        D_z_mean_face=field.D_z_mean_face,
        D_z_mean_cell=D_mean_cell,
        gauss_relative_residual=field.gauss_relative_residual,
        constitutive_flux_mismatch=mismatch,
        effective_step_eps_r=eps_step,
        memory_polarization=P_memory,
    )
