from __future__ import annotations

import math
import numpy as np

from pvdf_pf.morphology.three_phase import CRYSTAL, OAF, MAF


def _validate_fractions(crystal_fraction: float, oaf_fraction: float, maf_fraction: float) -> None:
    values = np.asarray([crystal_fraction, oaf_fraction, maf_fraction], dtype=float)
    if np.any(values <= 0.0):
        raise ValueError("all three phase fractions must be positive")
    if not np.isclose(values.sum(), 1.0, atol=1e-10):
        raise ValueError("phase fractions must sum to one")


def winding_normal_angle_deg(grid, winding_z: int, winding_x: int) -> float:
    """Return the physical normal angle measured from +z toward +x.

    The phase coordinate is periodic because the integer winding numbers multiply
    normalized cell coordinates.  This avoids a non-periodic seam at the cell edge.
    """
    if winding_z == 0 and winding_x == 0:
        raise ValueError("at least one winding number must be non-zero")
    kz = winding_z / (grid.nz * grid.dz)
    kx = winding_x / (grid.nx * grid.dx)
    return float(math.degrees(math.atan2(kx, kz)))


def periodic_winding_three_phase(
    grid,
    *,
    crystal_fraction: float,
    oaf_fraction: float,
    maf_fraction: float,
    winding_z: int = 1,
    winding_x: int = 0,
    offset: float = 0.0,
) -> np.ndarray:
    """Generate a periodic three-phase laminate with arbitrary commensurate tilt.

    The scalar phase coordinate is

        u = (m_z z/L_z + m_x x/L_x - offset) mod 1.

    Crystal occupies the first fraction of each period.  OAF is split equally on
    both sides of the crystal interval so that it represents an interfacial shell;
    MAF fills the remaining interval.

    This is a morphology-family generator, not a claim that PVDF lamellae are
    perfectly periodic or have a single orientation.
    """
    _validate_fractions(crystal_fraction, oaf_fraction, maf_fraction)
    if winding_z == 0 and winding_x == 0:
        raise ValueError("at least one winding number must be non-zero")

    iz = np.arange(grid.nz, dtype=float)[:, None] / grid.nz
    ix = np.arange(grid.nx, dtype=float)[None, :] / grid.nx
    u = (winding_z * iz + winding_x * ix - float(offset)) % 1.0

    half_oaf = 0.5 * oaf_fraction
    c0 = half_oaf
    c1 = c0 + crystal_fraction

    phase = np.full(grid.shape, MAF, dtype=np.int8)
    crystal = (u >= c0) & (u < c1)
    oaf_left = u < c0
    oaf_right = (u >= c1) & (u < c1 + half_oaf)
    phase[crystal] = CRYSTAL
    phase[oaf_left | oaf_right] = OAF
    return phase


def wavy_winding_three_phase(
    grid,
    *,
    crystal_fraction: float,
    oaf_fraction: float,
    maf_fraction: float,
    winding_z: int = 1,
    winding_x: int = 0,
    waviness_amplitude: float = 0.05,
    waviness_mode_z: int = 0,
    waviness_mode_x: int = 1,
    waviness_phase: float = 0.0,
    offset: float = 0.0,
) -> np.ndarray:
    """Generate a periodic sinusoidally distorted three-phase laminate.

    `waviness_amplitude` is expressed in units of one laminate period.  The
    perturbation is periodic and therefore compatible with the electrostatic cell
    problem.  It is used only to test morphology sensitivity in v0.1.
    """
    _validate_fractions(crystal_fraction, oaf_fraction, maf_fraction)
    if winding_z == 0 and winding_x == 0:
        raise ValueError("at least one winding number must be non-zero")
    if waviness_amplitude < 0.0:
        raise ValueError("waviness_amplitude must be non-negative")
    if waviness_mode_z == 0 and waviness_mode_x == 0:
        raise ValueError("waviness mode must contain at least one non-zero component")

    iz = np.arange(grid.nz, dtype=float)[:, None] / grid.nz
    ix = np.arange(grid.nx, dtype=float)[None, :] / grid.nx
    perturb = waviness_amplitude * np.sin(
        2.0 * np.pi * (waviness_mode_z * iz + waviness_mode_x * ix) + float(waviness_phase)
    )
    u = (winding_z * iz + winding_x * ix + perturb - float(offset)) % 1.0

    half_oaf = 0.5 * oaf_fraction
    c0 = half_oaf
    c1 = c0 + crystal_fraction

    phase = np.full(grid.shape, MAF, dtype=np.int8)
    crystal = (u >= c0) & (u < c1)
    oaf = (u < c0) | ((u >= c1) & (u < c1 + half_oaf))
    phase[crystal] = CRYSTAL
    phase[oaf] = OAF
    return phase
