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


def phase_fractions(phase_map: np.ndarray) -> dict[str, float]:
    n = phase_map.size
    return {
        "crystal": float(np.count_nonzero(phase_map == CRYSTAL) / n),
        "oaf": float(np.count_nonzero(phase_map == OAF) / n),
        "maf": float(np.count_nonzero(phase_map == MAF) / n),
    }
