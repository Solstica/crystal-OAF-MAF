import numpy as np

CRYSTAL = 0
OAF = 1
MAF = 2
PHASE_NAMES = {CRYSTAL: "crystal", OAF: "OAF", MAF: "MAF"}


def lamellar_three_phase(
    grid,
    crystal_thickness: int = 10,
    oaf_thickness: int = 3,
    period: int = 32,
    offset: int = 0,
) -> np.ndarray:
    """
    Deterministic periodic lamellar morphology for numerical verification.

    This geometry generator is not yet a statistical model of measured PVDF lamellae.
    """
    if crystal_thickness <= 0 or oaf_thickness <= 0:
        raise ValueError("crystal_thickness and oaf_thickness must be positive")
    if crystal_thickness + 2 * oaf_thickness >= period:
        raise ValueError("period must leave a non-zero MAF region")

    z = np.arange(grid.nz)[:, None]
    pos = (z - offset) % period
    phase = np.full(grid.shape, MAF, dtype=np.int8)

    crystal = pos < crystal_thickness
    oaf_top = (pos >= crystal_thickness) & (pos < crystal_thickness + oaf_thickness)
    oaf_bottom = pos >= (period - oaf_thickness)

    phase[np.broadcast_to(crystal, grid.shape)] = CRYSTAL
    phase[np.broadcast_to(oaf_top | oaf_bottom, grid.shape)] = OAF
    return phase


def lamellar_three_phase_from_fractions(
    grid,
    *,
    crystal_fraction: float,
    oaf_fraction: float,
    maf_fraction: float | None = None,
    period: int = 100,
    offset: int = 0,
) -> np.ndarray:
    """Construct a periodic crystal/OAF/MAF laminate from target volume fractions.

    OAF is split equally across the two crystal-amorphous interfaces of each
    period.  Fractions are converted to pixel occupancy, so the realized values
    converge to the targets as the period resolution increases.
    """
    fc = float(crystal_fraction)
    fo = float(oaf_fraction)
    fm = 1.0 - fc - fo if maf_fraction is None else float(maf_fraction)
    if min(fc, fo, fm) <= 0.0:
        raise ValueError("all three phase fractions must be positive")
    if not np.isclose(fc + fo + fm, 1.0, atol=1e-10):
        raise ValueError("phase fractions must sum to one")
    if period < 4:
        raise ValueError("period must be at least 4 pixels")

    z = np.arange(grid.nz)[:, None]
    u = ((z - offset) % period) / float(period)
    phase = np.full(grid.shape, MAF, dtype=np.int8)

    crystal = u < fc
    oaf_upper = (u >= fc) & (u < fc + 0.5 * fo)
    oaf_lower = u >= (1.0 - 0.5 * fo)
    phase[np.broadcast_to(crystal, grid.shape)] = CRYSTAL
    phase[np.broadcast_to(oaf_upper | oaf_lower, grid.shape)] = OAF
    return phase


def phase_fractions(phase_map: np.ndarray) -> dict[str, float]:
    n = phase_map.size
    return {
        "crystal": float(np.count_nonzero(phase_map == CRYSTAL) / n),
        "oaf": float(np.count_nonzero(phase_map == OAF) / n),
        "maf": float(np.count_nonzero(phase_map == MAF) / n),
    }
