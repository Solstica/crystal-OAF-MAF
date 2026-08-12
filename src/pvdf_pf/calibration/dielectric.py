"""Linear-dielectric cell problems for v0.1.

These utilities are deliberately separate from the TDGL `eps_b` map. A measured
small-signal relative permittivity contains dipolar response that may later be
represented explicitly by the polarization order parameter; using it directly as a
TDGL background permittivity would double count that response.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq
from scipy.sparse.linalg import LinearOperator, cg

from pvdf_pf.morphology.three_phase import CRYSTAL, OAF, MAF


def _harmonic_face(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Harmonic face value for flux continuity across sharp dielectric interfaces."""
    return 2.0 * a * b / (a + b)


def _neg_div_eps_grad(phi: np.ndarray, eps: np.ndarray, dz: float, dx: float) -> np.ndarray:
    eps_zp = _harmonic_face(eps, np.roll(eps, -1, axis=0))
    eps_zm = _harmonic_face(eps, np.roll(eps, 1, axis=0))
    eps_xp = _harmonic_face(eps, np.roll(eps, -1, axis=1))
    eps_xm = _harmonic_face(eps, np.roll(eps, 1, axis=1))

    gp_z = (np.roll(phi, -1, axis=0) - phi) / dz
    gm_z = (phi - np.roll(phi, 1, axis=0)) / dz
    gp_x = (np.roll(phi, -1, axis=1) - phi) / dx
    gm_x = (phi - np.roll(phi, 1, axis=1)) / dx

    return -(
        (eps_zp * gp_z - eps_zm * gm_z) / dz
        + (eps_xp * gp_x - eps_xm * gm_x) / dx
    )


def _macro_flux_divergence(eps: np.ndarray, grid, axis: str, E0: float) -> np.ndarray:
    if axis == "z":
        eps_p = _harmonic_face(eps, np.roll(eps, -1, axis=0))
        eps_m = _harmonic_face(eps, np.roll(eps, 1, axis=0))
        return E0 * (eps_p - eps_m) / grid.dz
    if axis == "x":
        eps_p = _harmonic_face(eps, np.roll(eps, -1, axis=1))
        eps_m = _harmonic_face(eps, np.roll(eps, 1, axis=1))
        return E0 * (eps_p - eps_m) / grid.dx
    raise ValueError("axis must be 'z' or 'x'")


def effective_permittivity(
    eps_r: np.ndarray,
    grid,
    *,
    axis: str = "z",
    E0: float = 1.0,
    tol: float = 1e-10,
    maxiter: int = 1000,
) -> tuple[float, np.ndarray]:
    r"""Return the periodic-cell effective relative permittivity.

    The corrector potential psi solves

        div[eps_r(r) (E0 e_axis - grad psi)] = 0.

    The returned effective permittivity is <D_axis>/E0 in relative units.
    """
    eps = np.asarray(eps_r, dtype=float)
    if eps.shape != grid.shape:
        raise ValueError("eps_r shape must match grid")
    if np.any(eps <= 0.0):
        raise ValueError("relative permittivity must be positive")
    if E0 == 0.0:
        raise ValueError("E0 must be non-zero")

    rhs = -_macro_flux_divergence(eps, grid, axis, E0)
    rhs -= rhs.mean()
    shape = grid.shape
    n = eps.size
    gauge = max(float(np.mean(eps)), 1.0) * 1e-12

    def matvec(x: np.ndarray) -> np.ndarray:
        psi = x.reshape(shape)
        y = _neg_div_eps_grad(psi, eps, grid.dz, grid.dx)
        return (y + gauge * psi.mean()).ravel()

    operator = LinearOperator((n, n), matvec=matvec, dtype=float)
    psi, info = cg(operator, rhs.ravel(), rtol=tol, atol=0.0, maxiter=maxiter)
    if info != 0:
        raise RuntimeError(f"dielectric cell-problem CG did not converge, info={info}")

    psi = psi.reshape(shape)
    psi -= psi.mean()

    if axis == "z":
        eps_face = _harmonic_face(eps, np.roll(eps, -1, axis=0))
        grad_face = (np.roll(psi, -1, axis=0) - psi) / grid.dz
    else:
        eps_face = _harmonic_face(eps, np.roll(eps, -1, axis=1))
        grad_face = (np.roll(psi, -1, axis=1) - psi) / grid.dx

    D_face = eps_face * (E0 - grad_face)
    eps_eff = float(D_face.mean() / E0)
    return eps_eff, psi


def phase_eps_map(phase_map: np.ndarray, phase_eps: dict[str, float]) -> np.ndarray:
    required = {"crystal", "oaf", "maf"}
    missing = required - set(phase_eps)
    if missing:
        raise KeyError(f"missing phase permittivities: {sorted(missing)}")
    out = np.empty(phase_map.shape, dtype=float)
    out[phase_map == CRYSTAL] = phase_eps["crystal"]
    out[phase_map == OAF] = phase_eps["oaf"]
    out[phase_map == MAF] = phase_eps["maf"]
    return out


def infer_unknown_phase_permittivity(
    phase_map: np.ndarray,
    grid,
    *,
    target_eps_eff: float,
    known_phase_eps: dict[str, float],
    unknown_phase: str = "oaf",
    axis: str = "z",
    bracket: tuple[float, float] = (1.01, 500.0),
) -> float:
    """Infer one effective local phase permittivity under an explicit morphology hypothesis.

    The result is an inverse-model parameter, not a direct experimental measurement.
    If the target cannot be reached within the bracket, a ValueError reports the
    attainable effective-permittivity interval.
    """
    if unknown_phase not in {"crystal", "oaf", "maf"}:
        raise ValueError("unknown_phase must be crystal, oaf, or maf")
    if unknown_phase in known_phase_eps:
        raise ValueError(f"{unknown_phase} must be omitted from known_phase_eps")

    lo, hi = map(float, bracket)
    if not (0.0 < lo < hi):
        raise ValueError("invalid positive bracket")

    def residual(value: float) -> float:
        eps_dict = dict(known_phase_eps)
        eps_dict[unknown_phase] = value
        eps_eff, _ = effective_permittivity(phase_eps_map(phase_map, eps_dict), grid, axis=axis)
        return eps_eff - target_eps_eff

    r_lo = residual(lo)
    r_hi = residual(hi)
    if r_lo == 0.0:
        return lo
    if r_hi == 0.0:
        return hi
    if r_lo * r_hi > 0.0:
        eff_lo = r_lo + target_eps_eff
        eff_hi = r_hi + target_eps_eff
        raise ValueError(
            f"target eps_eff={target_eps_eff:g} is outside the attainable interval "
            f"[{min(eff_lo, eff_hi):.6g}, {max(eff_lo, eff_hi):.6g}] for axis={axis!r} "
            f"and bracket={bracket}"
        )
    return float(brentq(residual, lo, hi, xtol=1e-10, rtol=1e-10, maxiter=200))
