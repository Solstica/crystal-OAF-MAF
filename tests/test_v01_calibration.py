import numpy as np

from pvdf_pf.calibration.dielectric import effective_permittivity
from pvdf_pf.calibration.polarization import oaf_polarization_lower_bound
from pvdf_pf.calibration.targets import validate_fraction_closure
from pvdf_pf.core.grid import Grid2D
from pvdf_pf.morphology.three_phase import lamellar_three_phase_from_fractions, phase_fractions


def test_fraction_closure_and_realized_morphology():
    assert validate_fraction_closure() == 1.0
    grid = Grid2D(nz=200, nx=8)
    phase = lamellar_three_phase_from_fractions(
        grid,
        crystal_fraction=0.52,
        oaf_fraction=0.28,
        maf_fraction=0.20,
        period=100,
    )
    fr = phase_fractions(phase)
    assert fr == {"crystal": 0.52, "oaf": 0.28, "maf": 0.20}


def test_uniform_dielectric_cell_problem():
    grid = Grid2D(nz=32, nx=24)
    eps = np.full(grid.shape, 7.5)
    for axis in ("x", "z"):
        value, _ = effective_permittivity(eps, grid, axis=axis)
        assert np.isclose(value, 7.5, rtol=1e-10, atol=1e-10)


def test_layered_dielectric_recovers_parallel_and_series_limits():
    grid = Grid2D(nz=100, nx=12)
    eps = np.full(grid.shape, 2.0)
    eps[50:, :] = 8.0

    eps_parallel, _ = effective_permittivity(eps, grid, axis="x")
    eps_series, _ = effective_permittivity(eps, grid, axis="z")

    arithmetic = 0.5 * 2.0 + 0.5 * 8.0
    harmonic = 1.0 / (0.5 / 2.0 + 0.5 / 8.0)
    assert np.isclose(eps_parallel, arithmetic, rtol=1e-8, atol=1e-8)
    assert np.isclose(eps_series, harmonic, rtol=2e-2, atol=2e-2)


def test_oaf_polarization_bound_is_physical():
    bound = oaf_polarization_lower_bound(
        film_ps_Cpm2=0.140,
        crystal_fraction=0.52,
        oaf_fraction=0.28,
        crystal_ps_upper_bound_Cpm2=0.188,
        maf_fraction=0.20,
        maf_ps_Cpm2=0.0,
    )
    assert np.isclose(bound, 0.15085714285714288)
