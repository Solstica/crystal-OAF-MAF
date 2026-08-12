import numpy as np

from pvdf_pf.core.grid import Grid2D
from pvdf_pf.morphology.oriented import (
    periodic_winding_three_phase,
    winding_normal_angle_deg,
    wavy_winding_three_phase,
)
from pvdf_pf.morphology.three_phase import phase_fractions


def test_periodic_winding_family_preserves_three_phases_and_fraction_scale():
    grid = Grid2D(nz=200, nx=200)
    phase = periodic_winding_three_phase(
        grid,
        crystal_fraction=0.52,
        oaf_fraction=0.28,
        maf_fraction=0.20,
        winding_z=1,
        winding_x=1,
    )
    frac = phase_fractions(phase)
    assert set(np.unique(phase)) == {0, 1, 2}
    assert abs(frac["crystal"] - 0.52) < 0.01
    assert abs(frac["oaf"] - 0.28) < 0.01
    assert abs(frac["maf"] - 0.20) < 0.01


def test_winding_angle_on_square_grid():
    grid = Grid2D(nz=100, nx=100, dz=1.0, dx=1.0)
    assert np.isclose(winding_normal_angle_deg(grid, 1, 0), 0.0)
    assert np.isclose(winding_normal_angle_deg(grid, 1, 1), 45.0)
    assert np.isclose(winding_normal_angle_deg(grid, 0, 1), 90.0)


def test_wavy_morphology_is_periodic_grid_compatible_and_three_phase():
    grid = Grid2D(nz=128, nx=128)
    phase = wavy_winding_three_phase(
        grid,
        crystal_fraction=0.52,
        oaf_fraction=0.28,
        maf_fraction=0.20,
        winding_z=1,
        winding_x=0,
        waviness_amplitude=0.08,
        waviness_mode_x=2,
    )
    assert phase.shape == grid.shape
    assert set(np.unique(phase)) == {0, 1, 2}
