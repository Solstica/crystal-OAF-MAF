import numpy as np


def grad_periodic(f: np.ndarray, dz: float, dx: float) -> tuple[np.ndarray, np.ndarray]:
    """Second-order centered periodic gradient, returned as (df/dz, df/dx)."""
    df_dz = (np.roll(f, -1, axis=0) - np.roll(f, 1, axis=0)) / (2.0 * dz)
    df_dx = (np.roll(f, -1, axis=1) - np.roll(f, 1, axis=1)) / (2.0 * dx)
    return df_dz, df_dx


def semi_implicit_gradient_step(
    P: np.ndarray,
    dF_nongrad: np.ndarray,
    grid,
    mobility: float,
    dt: float,
    kappa: float,
) -> np.ndarray:
    """One scalar TDGL step with constant isotropic gradient energy treated implicitly."""
    _, _, k2 = grid.kgrid()
    rhs = P - dt * mobility * dF_nongrad
    rhs_k = np.fft.rfftn(rhs)
    denom = 1.0 + dt * mobility * kappa * k2
    P_next = np.fft.irfftn(rhs_k / denom, s=grid.shape)
    return np.asarray(P_next.real, dtype=float)
