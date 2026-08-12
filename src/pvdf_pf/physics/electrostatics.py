import numpy as np
from scipy.sparse.linalg import LinearOperator, cg

from pvdf_pf.core.spectral import grad_periodic


def _harmonic_face(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Harmonic face value for normal dielectric flux across sharp interfaces."""
    return 2.0 * a * b / (a + b)


def _neg_div_eps_grad(phi: np.ndarray, eps: np.ndarray, dz: float, dx: float) -> np.ndarray:
    """Periodic finite-volume operator -div(eps grad(phi)) using harmonic face eps."""
    eps_zp = _harmonic_face(eps, np.roll(eps, -1, axis=0))
    eps_zm = _harmonic_face(eps, np.roll(eps, 1, axis=0))
    gp_z = (np.roll(phi, -1, axis=0) - phi) / dz
    gm_z = (phi - np.roll(phi, 1, axis=0)) / dz
    term_z = -(eps_zp * gp_z - eps_zm * gm_z) / dz

    eps_xp = _harmonic_face(eps, np.roll(eps, -1, axis=1))
    eps_xm = _harmonic_face(eps, np.roll(eps, 1, axis=1))
    gp_x = (np.roll(phi, -1, axis=1) - phi) / dx
    gm_x = (phi - np.roll(phi, 1, axis=1)) / dx
    term_x = -(eps_xp * gp_x - eps_xm * gm_x) / dx
    return term_z + term_x


def solve_depolarization_scalar(
    Pz: np.ndarray,
    eps_b: np.ndarray,
    grid,
    *,
    eps0: float = 1.0,
    tol: float = 1e-9,
    maxiter: int = 500,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r"""
    Solve the periodic heterogeneous electrostatic equation

        div(eps0 * eps_b(r) * grad(phi)) = dPz/dz,
        E_dep = -grad(phi).

    A tiny mean-projector term removes the constant-potential null mode without
    changing non-zero Fourier components.
    """
    eps = eps0 * np.asarray(eps_b, dtype=float)
    if eps.shape != grid.shape:
        raise ValueError("eps_b shape must match grid")
    if np.any(eps <= 0.0):
        raise ValueError("eps0 * eps_b must be positive")

    dP_dz = (np.roll(Pz, -1, axis=0) - np.roll(Pz, 1, axis=0)) / (2.0 * grid.dz)
    rhs = -dP_dz
    rhs -= rhs.mean()

    shape = grid.shape
    n = Pz.size
    gauge = max(float(np.mean(eps)), 1.0) * 1e-12

    def matvec(x: np.ndarray) -> np.ndarray:
        phi = x.reshape(shape)
        y = _neg_div_eps_grad(phi, eps, grid.dz, grid.dx)
        y = y + gauge * phi.mean()
        return y.ravel()

    operator = LinearOperator((n, n), matvec=matvec, dtype=float)
    phi, info = cg(operator, rhs.ravel(), rtol=tol, atol=0.0, maxiter=maxiter)
    if info != 0:
        raise RuntimeError(f"heterogeneous Poisson CG did not converge, info={info}")

    phi = phi.reshape(shape)
    phi -= phi.mean()
    dphi_dz, dphi_dx = grad_periodic(phi, grid.dz, grid.dx)
    return phi, -dphi_dz, -dphi_dx
