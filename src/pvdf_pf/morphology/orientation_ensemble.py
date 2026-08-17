"""Same-RMS orientation-spectrum morphology utilities for v0.1.19.

The v0.1.18 descriptor ``sqrt(<n_ND^2>)`` constrains only the second moment of
local interface orientation.  This module constructs several smooth periodic
interface spectra with the same RMS normal projection so the electrostatic model
can test whether that one descriptor is sufficient.

The harmonic spectra defined here are project geometry hypotheses.  Their mode
numbers, relative amplitudes and phases are not inferred from Huang et al.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from pvdf_pf.morphology.three_phase import CRYSTAL, OAF, MAF


@dataclass(frozen=True)
class HarmonicMode:
    mode_z: int
    coefficient: float
    phase_rad: float = 0.0

    def validate(self) -> None:
        if self.mode_z <= 0:
            raise ValueError("mode_z must be a positive integer")
        if not np.isfinite(self.coefficient) or self.coefficient == 0.0:
            raise ValueError("coefficient must be finite and non-zero")
        if not np.isfinite(self.phase_rad):
            raise ValueError("phase_rad must be finite")


def _validate_fractions(crystal_fraction: float, oaf_fraction: float, iaf_fraction: float) -> None:
    values = np.asarray([crystal_fraction, oaf_fraction, iaf_fraction], dtype=float)
    if np.any(values <= 0.0):
        raise ValueError("all structural fractions must be positive")
    if not np.isclose(values.sum(), 1.0, atol=1e-10):
        raise ValueError("structural fractions must sum to one")


def spectral_phase_coordinate_and_normal(
    grid,
    modes: tuple[HarmonicMode, ...] | list[HarmonicMode],
    scale: float,
    *,
    offset: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Return periodic phase coordinate ``u`` and its film-normal normal component.

    The mean laminate coordinate is ``x/Lx``.  A z-dependent smooth periodic
    displacement is added,

        u = x/Lx + scale * sum_k c_k sin(2*pi*k*z/Lz + phi_k).

    The local interface normal is parallel to ``grad(u)``.  Because the
    perturbation depends only on z, ``du/dx = 1/Lx`` remains positive and every
    interface stays a single-valued graph in x.
    """
    if not np.isfinite(scale) or scale < 0.0:
        raise ValueError("scale must be finite and non-negative")
    modes = tuple(modes)
    if not modes:
        raise ValueError("at least one harmonic mode is required")
    for mode in modes:
        mode.validate()

    zeta = np.arange(grid.nz, dtype=float)[:, None] / grid.nz
    xi = np.arange(grid.nx, dtype=float)[None, :] / grid.nx
    perturb = np.zeros((grid.nz, 1), dtype=float)
    dperturb_dz = np.zeros((grid.nz, 1), dtype=float)
    Lz = float(grid.nz * grid.dz)
    Lx = float(grid.nx * grid.dx)

    for mode in modes:
        arg = 2.0 * np.pi * mode.mode_z * zeta + float(mode.phase_rad)
        perturb += float(mode.coefficient) * np.sin(arg)
        dperturb_dz += (
            float(mode.coefficient)
            * 2.0
            * np.pi
            * mode.mode_z
            * np.cos(arg)
            / Lz
        )

    perturb *= float(scale)
    dperturb_dz *= float(scale)
    u = (xi + perturb - float(offset)) % 1.0

    du_dx = 1.0 / Lx
    denom = np.sqrt(dperturb_dz**2 + du_dx**2)
    n_nd_column = dperturb_dz / denom
    n_nd = np.broadcast_to(n_nd_column, grid.shape).copy()
    return u, n_nd


def spectral_orientation_metrics(
    grid,
    modes: tuple[HarmonicMode, ...] | list[HarmonicMode],
    scale: float,
) -> dict[str, float]:
    _, n_nd = spectral_phase_coordinate_and_normal(grid, modes, scale)
    abs_n = np.abs(n_nd)
    return {
        "rms_n_ND": float(np.sqrt(np.mean(n_nd**2))),
        "mean_abs_n_ND": float(np.mean(abs_n)),
        "max_abs_n_ND": float(np.max(abs_n)),
        "mean_n_ND": float(np.mean(n_nd)),
        "n_ND_p95_abs": float(np.quantile(abs_n, 0.95)),
        "n_ND_p99_abs": float(np.quantile(abs_n, 0.99)),
    }


def scale_for_target_rms_normal(
    grid,
    modes: tuple[HarmonicMode, ...] | list[HarmonicMode],
    target_rms_n_ND: float,
    *,
    scale_bracket: tuple[float, float] = (0.0, 1.0),
) -> float:
    """Scale one harmonic spectrum to a requested RMS film-normal projection."""
    target = float(target_rms_n_ND)
    if not (0.0 <= target < 1.0):
        raise ValueError("target_rms_n_ND must lie in [0,1)")
    if target == 0.0:
        return 0.0

    lo, hi = map(float, scale_bracket)
    if lo < 0.0 or hi <= lo:
        raise ValueError("invalid scale bracket")

    def residual(scale: float) -> float:
        return spectral_orientation_metrics(grid, modes, scale)["rms_n_ND"] - target

    r_lo = residual(lo)
    r_hi = residual(hi)
    if r_lo > 0.0 or r_hi < 0.0:
        attained = (
            spectral_orientation_metrics(grid, modes, lo)["rms_n_ND"],
            spectral_orientation_metrics(grid, modes, hi)["rms_n_ND"],
        )
        raise ValueError(
            f"target_rms_n_ND={target:g} outside attainable interval "
            f"[{attained[0]:.6g}, {attained[1]:.6g}]"
        )
    return float(brentq(residual, lo, hi, xtol=1e-12, rtol=1e-12, maxiter=200))


def spectral_three_phase(
    grid,
    *,
    crystal_fraction: float,
    oaf_fraction: float,
    iaf_fraction: float,
    modes: tuple[HarmonicMode, ...] | list[HarmonicMode],
    scale: float,
    offset: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Build beta/OAF/IAF labels and return the local ``n_ND`` field."""
    _validate_fractions(crystal_fraction, oaf_fraction, iaf_fraction)
    u, n_nd = spectral_phase_coordinate_and_normal(
        grid, modes, scale, offset=offset
    )

    half_oaf = 0.5 * float(oaf_fraction)
    c0 = half_oaf
    c1 = c0 + float(crystal_fraction)
    phase = np.full(grid.shape, MAF, dtype=np.int8)
    crystal = (u >= c0) & (u < c1)
    oaf = (u < c0) | ((u >= c1) & (u < c1 + half_oaf))
    phase[crystal] = CRYSTAL
    phase[oaf] = OAF
    return phase, n_nd
