import numpy as np

from pvdf_pf.core.grid import Grid2D
from pvdf_pf.physics.electrostatics import solve_depolarization_scalar


def test_uniform_polarization_has_zero_periodic_depolarization_field():
    grid = Grid2D(24, 20)
    P = np.ones(grid.shape)
    eps = np.full(grid.shape, 4.0)
    phi, Ez, Ex = solve_depolarization_scalar(P, eps, grid)
    assert np.max(np.abs(phi)) < 1e-10
    assert np.max(np.abs(Ez)) < 1e-10
    assert np.max(np.abs(Ex)) < 1e-10
