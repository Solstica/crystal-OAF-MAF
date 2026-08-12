from pvdf_pf.core.grid import Grid2D
from pvdf_pf.morphology.three_phase import lamellar_three_phase, phase_fractions


def test_three_phase_present():
    grid = Grid2D(48, 32)
    phase = lamellar_three_phase(grid, crystal_thickness=8, oaf_thickness=2, period=16)
    assert set(phase.ravel()) == {0, 1, 2}
    fractions = phase_fractions(phase)
    assert abs(sum(fractions.values()) - 1.0) < 1e-12
