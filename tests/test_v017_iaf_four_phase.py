from pathlib import Path

import numpy as np

from pvdf_pf.calibration.four_phase_dielectric import build_four_phase_permittivity_state, solve_four_phase_cell
from pvdf_pf.calibration.iaf_permittivity import epsilon_iaf_from_melt_extrapolation, fit_molten_pvdf_permittivity
from pvdf_pf.calibration.project_devitrification import build_regularized_oaf_state
from pvdf_pf.core.grid import Grid2D
from pvdf_pf.morphology.mobility_scores import iaf_proximal_score
from pvdf_pf.morphology.oaf_devitrification import allocate_source_oaf_subphases, source_phase_fractions
from pvdf_pf.morphology.oriented import periodic_winding_three_phase

ROOT = Path(__file__).resolve().parents[1]
RUI = ROOT / "data/literature/rui2022"
FIG1B = RUI / "figure_1b_melt_eps_digitized.csv"
FIG5B = RUI / "figure_5b_maintext_digitized.csv"


def test_figure1b_linear_fit_and_low_temperature_extrapolation():
    fit = fit_molten_pvdf_permittivity(FIG1B)
    assert fit.n_points == 11
    assert fit.r_squared > 0.995
    assert -0.045 < fit.slope_per_C < -0.041
    assert 17.0 < fit.intercept < 18.0

    eps_m30, _ = epsilon_iaf_from_melt_extrapolation(FIG1B, temperature_C=-30.0)
    eps_40, _ = epsilon_iaf_from_melt_extrapolation(FIG1B, temperature_C=40.0)
    assert 18.0 < eps_m30 < 19.5
    assert 15.0 < eps_40 < 16.5
    assert eps_m30 > eps_40


def test_four_phase_parallel_cell_closes_and_is_uniform_field():
    grid = Grid2D(nz=50, nx=50)
    phase = periodic_winding_three_phase(
        grid,
        crystal_fraction=0.60,
        oaf_fraction=0.20,
        maf_fraction=0.20,
        winding_z=1,
        winding_x=0,
    )
    dev = build_regularized_oaf_state(FIG5B, temperature_C=20.0)
    state = dev.as_source_fraction_state(sample_state="melt-recrystallized")
    source_map, _ = allocate_source_oaf_subphases(
        phase,
        fraction_state=state,
        mobility_score=iaf_proximal_score(phase),
    )
    eps_state = build_four_phase_permittivity_state(
        RUI,
        temperature_C=20.0,
        sample_state="poled",
    )
    _, _, _, summary = solve_four_phase_cell(source_map, eps_state, grid, axis="x", E0=1.0)
    fractions = source_phase_fractions(source_map)

    expected = (
        fractions["crystal"] * eps_state.epsilon_crystal
        + fractions["roaf"] * eps_state.epsilon_roaf
        + fractions["moaf"] * eps_state.epsilon_moaf
        + fractions["iaf"] * eps_state.epsilon_iaf
    )
    assert np.isclose(summary.epsilon_effective, expected, rtol=2e-3, atol=2e-3)
    assert summary.field_std < 1e-8
    assert np.isclose(sum(fractions.values()), 1.0, atol=1e-12)
