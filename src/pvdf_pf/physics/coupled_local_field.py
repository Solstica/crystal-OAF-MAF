"""Self-consistent local-field coupling for generalized-Debye polarization.

v0.1.15 couples the physical-time auxiliary polarization bank to Gauss' law while
keeping the ferroelectric TDGL order parameter frozen.  The source response is still
the same-state combined OAF+IAF amorphous response established in v0.1.13-v0.1.14.

For one zero-order-hold step, each Debye mode obeys

    P_j^n = a_j P_j^(n-1) + (1-a_j) eps0 Delta_eps_j E^n,
    a_j = exp(-dt/tau_j).

Substituting this relation into Gauss' law gives one *linear* heterogeneous Poisson
problem per time step,

    div{ eps0 [eps_b + sum_j (1-a_j) Delta_eps_j] E^n
         + sum_j a_j P_j^(n-1) } = 0.

Thus the local field and the new auxiliary polarizations are self-consistent without
a fixed-point iteration.  This is a numerical constitutive coupling; it does not
calibrate the dimensionless TDGL clock.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse.linalg import LinearOperator, cg

from pvdf_pf.core.spectral import grad_periodic
from pvdf_pf.physics.dipolar import EPS0
from pvdf_pf.physics.electrostatics import _harmonic_face, _neg_div_eps_grad
from pvdf_pf.physics.generalized_debye import GeneralizedDebyeBank


@dataclass(frozen=True)
class LocalFieldState:
    potential: np.ndarray
    E_z: np.ndarray
    E_x: np.ndarray
    D_z_mean_face: float
    gauss_relative_residual: float


@dataclass(frozen=True)
class CoupledRelaxationStep:
    P_modes: np.ndarray
    P_total: np.ndarray
    potential: np.ndarray
    E_z: np.ndarray
    E_x: np.ndarray
    D_z_mean_face: float
    D_z_mean_cell: float
    gauss_relative_residual: float
    constitutive_flux_mismatch: float
    effective_step_eps_r: np.ndarray
    memory_polarization: np.ndarray


def _divergence_of_uniform_macro_flux(
    eps_abs: np.ndarray,
    *,
    E_external_z: float,
    dz: float,
) -> np.ndarray:
    eps_zp = _harmonic_face(eps_abs, np.roll(eps_abs, -1, axis=0))
    eps_zm = _harmonic_face(eps_abs, np.roll(eps_abs, 1, axis=0))
    return float(E_external_z) * (eps_zp - eps_zm) / float(dz)


def _divergence_of_cell_pz(Pz: np.ndarray, *, dz: float) -> np.ndarray:
    """Finite-volume divergence using arithmetic face interpolation of Pz."""
    p_zp = 0.5 * (Pz + np.roll(Pz, -1, axis=0))
    p_zm = 0.5 * (Pz + np.roll(Pz, 1, axis=0))
    return (p_zp - p_zm) / float(dz)


def solve_periodic_local_field_scalar(
    Pz_explicit: np.ndarray,
    eps_r: np.ndarray,
    grid,
    *,
    E_external_z: float,
    eps0: float = EPS0,
    tol: float = 1e-10,
    maxiter: int = 4000,
) -> LocalFieldState:
    r"""Solve periodic Gauss electrostatics with heterogeneous eps and explicit Pz.

    The electric field is

        E = E_external_z e_z - grad(psi)

    and the discretized equation is

        div[eps0 eps_r E + Pz e_z] = 0.

    Unlike the older depolarization-only helper, this function includes the
    dielectric correction driven by a spatially varying eps_r under the applied
    macroscopic field.
    """
    eps_rel = np.asarray(eps_r, dtype=float)
    Pz = np.asarray(Pz_explicit, dtype=float)
    if eps_rel.shape != grid.shape or Pz.shape != grid.shape:
        raise ValueError("eps_r and Pz_explicit must match grid.shape")
    if np.any(~np.isfinite(eps_rel)) or np.any(eps_rel <= 0.0):
        raise ValueError("eps_r must be finite and positive")
    if np.any(~np.isfinite(Pz)):
        raise ValueError("Pz_explicit must be finite")
    if not np.isfinite(E_external_z):
        raise ValueError("E_external_z must be finite")
    if eps0 <= 0.0 or tol <= 0.0 or maxiter <= 0:
        raise ValueError("eps0, tol and maxiter must be positive")

    eps_abs = float(eps0) * eps_rel
    macro_div = _divergence_of_uniform_macro_flux(
        eps_abs,
        E_external_z=float(E_external_z),
        dz=grid.dz,
    )
    p_div = _divergence_of_cell_pz(Pz, dz=grid.dz)
    rhs = -macro_div - p_div
    rhs -= rhs.mean()

    shape = grid.shape
    n = Pz.size
    gauge = max(float(np.mean(eps_abs)), float(eps0)) * 1e-12

    def matvec(x: np.ndarray) -> np.ndarray:
        psi = x.reshape(shape)
        y = _neg_div_eps_grad(psi, eps_abs, grid.dz, grid.dx)
        return (y + gauge * psi.mean()).ravel()

    operator = LinearOperator((n, n), matvec=matvec, dtype=float)

    eps_zp = _harmonic_face(eps_abs, np.roll(eps_abs, -1, axis=0))
    eps_zm = _harmonic_face(eps_abs, np.roll(eps_abs, 1, axis=0))
    eps_xp = _harmonic_face(eps_abs, np.roll(eps_abs, -1, axis=1))
    eps_xm = _harmonic_face(eps_abs, np.roll(eps_abs, 1, axis=1))
    diag = (
        (eps_zp + eps_zm) / (grid.dz * grid.dz)
        + (eps_xp + eps_xm) / (grid.dx * grid.dx)
        + gauge / n
    )
    if np.any(~np.isfinite(diag)) or np.any(diag <= 0.0):
        raise RuntimeError("invalid local-field Jacobi diagonal")
    inv_diag = 1.0 / diag.ravel()
    preconditioner = LinearOperator(
        (n, n), matvec=lambda x: inv_diag * x, dtype=float
    )

    psi, info = cg(
        operator,
        rhs.ravel(),
        M=preconditioner,
        rtol=tol,
        atol=0.0,
        maxiter=maxiter,
    )
    if info != 0:
        residual = operator.matvec(psi) - rhs.ravel()
        rhs_norm = max(float(np.linalg.norm(rhs.ravel())), 1e-30)
        rel = float(np.linalg.norm(residual) / rhs_norm)
        raise RuntimeError(
            "self-consistent local-field CG did not converge, "
            f"info={info}, relative_residual={rel:.3e}"
        )

    psi = psi.reshape(shape)
    psi -= psi.mean()
    dpsi_dz, dpsi_dx = grad_periodic(psi, grid.dz, grid.dx)
    E_z = float(E_external_z) - dpsi_dz
    E_x = -dpsi_dx

    grad_z_face = (np.roll(psi, -1, axis=0) - psi) / grid.dz
    grad_x_face = (np.roll(psi, -1, axis=1) - psi) / grid.dx
    P_z_face = 0.5 * (Pz + np.roll(Pz, -1, axis=0))
    D_z_face = eps_zp * (float(E_external_z) - grad_z_face) + P_z_face
    D_x_face = eps_xp * (-grad_x_face)

    div_D = (
        (D_z_face - np.roll(D_z_face, 1, axis=0)) / grid.dz
        + (D_x_face - np.roll(D_x_face, 1, axis=1)) / grid.dx
    )
    flux_scale = max(
        float(np.sqrt(np.mean(D_z_face**2))) / max(float(grid.dz), 1e-30),
        float(np.sqrt(np.mean(D_x_face**2))) / max(float(grid.dx), 1e-30),
        1e-30,
    )
    gauss_rel = float(np.sqrt(np.mean(div_D**2)) / flux_scale)

    return LocalFieldState(
        potential=psi,
        E_z=E_z,
        E_x=E_x,
        D_z_mean_face=float(D_z_face.mean()),
        gauss_relative_residual=gauss_rel,
    )


def build_combined_amorphous_background_map(
    phase_map: np.ndarray,
    bank: GeneralizedDebyeBank,
    *,
    crystal_phase_id: int = 0,
    epsilon_crystal_background: float = 3.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Return eps_b map and combined-amorphous mask for v0.1.15.

    All non-crystal cells share the same epsilon_infinity and relaxation bank.
    Keeping OAF/MAF labels in the morphology does not imply that v0.1.15 has
    independently measured OAF and IAF constitutive laws.
    """
    bank.validate()
    phase = np.asarray(phase_map)
    if phase.ndim != 2:
        raise ValueError("phase_map must be two-dimensional")
    if epsilon_crystal_background <= 0.0:
        raise ValueError("epsilon_crystal_background must be positive")
    amorphous = phase != int(crystal_phase_id)
    if not np.any(amorphous) or not np.any(~amorphous):
        raise ValueError("phase_map must contain both crystal and amorphous cells")
    eps_b = np.full(phase.shape, float(bank.epsilon_infinity), dtype=float)
    eps_b[~amorphous] = float(epsilon_crystal_background)
    return eps_b, amorphous


def advance_self_consistent_generalized_debye(
    P_modes_old: np.ndarray,
    eps_background_r: np.ndarray,
    amorphous_mask: np.ndarray,
    grid,
    *,
    E_external_z: float,
    dt_s: float,
    bank: GeneralizedDebyeBank,
    eps0: float = EPS0,
    tol: float = 1e-10,
    maxiter: int = 4000,
) -> CoupledRelaxationStep:
    """Advance one self-consistent zero-order-hold generalized-Debye step."""
    bank.validate()
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive")
    eps_b = np.asarray(eps_background_r, dtype=float)
    mask = np.asarray(amorphous_mask, dtype=bool)
    if eps_b.shape != grid.shape or mask.shape != grid.shape:
        raise ValueError("eps_background_r and amorphous_mask must match grid.shape")
    if np.any(eps_b <= 0.0) or np.any(~np.isfinite(eps_b)):
        raise ValueError("eps_background_r must be finite and positive")

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
        P_memory,
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

    D_cell = float(eps0) * eps_b * field.E_z + P_total
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


def aligned_laminate_effective_permittivity(
    epsilon_crystal: complex | float,
    epsilon_amorphous: complex | float,
    *,
    crystal_fraction: float,
    field_parallel_to_layers: bool,
) -> complex:
    """Analytic two-phase aligned-laminate reference for geometry audits."""
    fc = float(crystal_fraction)
    if not (0.0 < fc < 1.0):
        raise ValueError("crystal_fraction must lie in (0,1)")
    ec = complex(epsilon_crystal)
    ea = complex(epsilon_amorphous)
    if ec.real <= 0.0 or ea.real <= 0.0:
        raise ValueError("phase permittivities must have positive real part")
    fa = 1.0 - fc
    if field_parallel_to_layers:
        return fc * ec + fa * ea
    return 1.0 / (fc / ec + fa / ea)
