"""Orientation-spread utilities for the v0.1.18 Huang-SAXS constraint.

The local interface normal for the v0.1.17 morphology

    u = x/Lx + A sin(2*pi*mz*z/Lz)

has film-normal component n_z.  v0.1.18 parameterizes the morphology by the
observable-like geometric descriptor sqrt(<n_z^2>) rather than by the internal
sinusoidal amplitude A.

The target interval is not a direct Huang et al. measurement.  It is an
image-derived project constraint recorded separately in
``data/literature/huang2021/README.md``.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq


def sinusoidal_normal_projection_metrics(
    grid,
    amplitude: float,
    *,
    mode_z: int = 1,
    mode_x: int = 0,
) -> dict[str, float]:
    """Return normal-projection metrics for the sinusoidal laminate coordinate."""
    if amplitude < 0.0:
        raise ValueError("amplitude must be non-negative")
    if mode_z == 0 and mode_x == 0:
        raise ValueError("at least one waviness mode must be non-zero")

    z = np.arange(grid.nz, dtype=float)[:, None] / grid.nz
    x = np.arange(grid.nx, dtype=float)[None, :] / grid.nx
    phase = 2.0 * np.pi * (mode_z * z + mode_x * x)
    Lz = float(grid.nz * grid.dz)
    Lx = float(grid.nx * grid.dx)

    du_dz = float(amplitude) * 2.0 * np.pi * mode_z * np.cos(phase) / Lz
    du_dx = 1.0 / Lx + float(amplitude) * 2.0 * np.pi * mode_x * np.cos(phase) / Lx
    denom = np.sqrt(du_dz**2 + du_dx**2)
    nz = np.divide(du_dz, denom, out=np.zeros_like(denom), where=denom > 0.0)

    return {
        "rms_nz": float(np.sqrt(np.mean(nz**2))),
        "mean_abs_nz": float(np.mean(np.abs(nz))),
        "max_abs_nz": float(np.max(np.abs(nz))),
    }


def amplitude_for_target_rms_nz(
    grid,
    target_rms_nz: float,
    *,
    mode_z: int = 1,
    mode_x: int = 0,
    amplitude_bracket: tuple[float, float] = (0.0, 0.5),
) -> float:
    """Invert the sinusoidal geometry so ``sqrt(<n_z^2>)`` is the control variable."""
    target = float(target_rms_nz)
    if not (0.0 <= target < 1.0):
        raise ValueError("target_rms_nz must lie in [0, 1)")
    if target == 0.0:
        return 0.0

    lo, hi = map(float, amplitude_bracket)
    if lo < 0.0 or hi <= lo:
        raise ValueError("invalid amplitude bracket")

    def residual(amplitude: float) -> float:
        return sinusoidal_normal_projection_metrics(
            grid, amplitude, mode_z=mode_z, mode_x=mode_x
        )["rms_nz"] - target

    r_lo = residual(lo)
    r_hi = residual(hi)
    if r_lo > 0.0 or r_hi < 0.0:
        attainable = (
            sinusoidal_normal_projection_metrics(grid, lo, mode_z=mode_z, mode_x=mode_x)["rms_nz"],
            sinusoidal_normal_projection_metrics(grid, hi, mode_z=mode_z, mode_x=mode_x)["rms_nz"],
        )
        raise ValueError(
            f"target_rms_nz={target:g} is outside attainable interval "
            f"[{attainable[0]:.6g}, {attainable[1]:.6g}]"
        )

    return float(brentq(residual, lo, hi, xtol=1e-12, rtol=1e-12, maxiter=200))
