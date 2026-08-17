import numpy as np

from pvdf_pf.calibration.dielectric import effective_permittivity
from pvdf_pf.core.grid import Grid2D
from pvdf_pf.morphology.oriented import periodic_winding_three_phase
from pvdf_pf.morphology.three_phase import CRYSTAL, phase_fractions
from pvdf_pf.physics.coupled_local_field import (
    advance_self_consistent_generalized_debye,
    aligned_laminate_effective_permittivity,
    build_combined_amorphous_background_map,
    solve_periodic_local_field_scalar,
)
from pvdf_pf.physics.dipolar import EPS0
from pvdf_pf.physics.generalized_debye import GeneralizedDebyeBank


def test_uniform_local_field_with_uniform_explicit_polarization():
    grid = Grid2D(nz=12, nx=10, dz=1.0, dx=1.0)
    eps_r = np.full(grid.shape, 7.0)
    Pz = np.full(grid.shape, 0.02)
    E0 = 2.5e5
    state = solve_periodic_local_field_scalar(Pz, eps_r, grid, E_external_z=E0)
    assert np.max(np.abs(state.potential)) < 1e-12
    assert np.allclose(state.E_z, E0, rtol=0, atol=1e-10)
    assert np.allclose(state.E_x, 0.0, rtol=0, atol=1e-12)
    assert np.isclose(state.D_z_mean_face, EPS0 * 7.0 * E0 + 0.02, rtol=1e-12)
    assert state.gauss_relative_residual < 1e-12


def test_uniform_coupled_step_matches_exact_zoh_relation():
    grid = Grid2D(nz=8, nx=9)
    bank = GeneralizedDebyeBank(
        epsilon_infinity=5.0,
        tau_modes_s=(0.2, 2.0),
        delta_epsilon_modes=(2.0, 4.0),
    )
    mask = np.ones(grid.shape, dtype=bool)
    eps_b = np.full(grid.shape, bank.epsilon_infinity)
    P0 = np.zeros((bank.n_modes,) + grid.shape)
    E0 = 1.2e6
    dt = 0.15

    step = advance_self_consistent_generalized_debye(
        P0,
        eps_b,
        mask,
        grid,
        E_external_z=E0,
        dt_s=dt,
        bank=bank,
    )
    tau = np.asarray(bank.tau_modes_s)
    strength = np.asarray(bank.delta_epsilon_modes)
    expected_modes = (1.0 - np.exp(-dt / tau)) * EPS0 * strength * E0
    assert np.allclose(step.E_z, E0, rtol=0, atol=1e-9)
    for j, expected in enumerate(expected_modes):
        assert np.allclose(step.P_modes[j], expected, rtol=1e-12, atol=1e-20)
    expected_eps_step = bank.epsilon_infinity + np.sum(
        (1.0 - np.exp(-dt / tau)) * strength
    )
    expected_D = EPS0 * expected_eps_step * E0
    assert np.isclose(step.D_z_mean_face, expected_D, rtol=1e-12)
    assert np.isclose(step.D_z_mean_cell, expected_D, rtol=1e-12)
    assert step.gauss_relative_residual < 1e-12
    assert step.constitutive_flux_mismatch < 1e-12


def test_quasistatic_coupled_limit_matches_static_cell_problem_for_both_orientations():
    grid = Grid2D(nz=50, nx=50)
    bank = GeneralizedDebyeBank(
        epsilon_infinity=6.0,
        tau_modes_s=(0.01, 0.1, 1.0),
        delta_epsilon_modes=(4.0, 8.0, 12.0),
    )
    fractions = dict(crystal_fraction=0.52, oaf_fraction=0.28, maf_fraction=0.20)
    E0 = 1.0e6
    dt = 100.0 * max(bank.tau_modes_s)

    for winding_z, winding_x, parallel in ((0, 1, True), (1, 0, False)):
        phase = periodic_winding_three_phase(
            grid, **fractions, winding_z=winding_z, winding_x=winding_x
        )
        eps_b, mask = build_combined_amorphous_background_map(
            phase, bank, epsilon_crystal_background=3.0
        )
        P0 = np.zeros((bank.n_modes,) + grid.shape)
        step = advance_self_consistent_generalized_debye(
            P0,
            eps_b,
            mask,
            grid,
            E_external_z=E0,
            dt_s=dt,
            bank=bank,
        )

        eps_static = np.full(grid.shape, bank.epsilon_infinity + bank.delta_epsilon_total)
        eps_static[phase == CRYSTAL] = 3.0
        eps_ref, _ = effective_permittivity(eps_static, grid, axis="z")
        eps_coupled = step.D_z_mean_face / (EPS0 * E0)
        eps_analytic = aligned_laminate_effective_permittivity(
            3.0,
            bank.epsilon_infinity + bank.delta_epsilon_total,
            crystal_fraction=0.52,
            field_parallel_to_layers=parallel,
        ).real

        assert np.isclose(eps_coupled, eps_ref, rtol=2e-8, atol=2e-8)
        assert np.isclose(eps_coupled, eps_analytic, rtol=2e-8, atol=2e-8)
        assert step.gauss_relative_residual < 1e-7


def test_project_phase_labels_realize_combined_amorphous_fraction_without_oaf_promotion():
    grid = Grid2D(nz=50, nx=50)
    phase = periodic_winding_three_phase(
        grid,
        crystal_fraction=0.52,
        oaf_fraction=0.28,
        maf_fraction=0.20,
        winding_z=0,
        winding_x=1,
    )
    bank = GeneralizedDebyeBank(
        epsilon_infinity=8.0,
        tau_modes_s=(1.0,),
        delta_epsilon_modes=(2.0,),
    )
    _, mask = build_combined_amorphous_background_map(phase, bank)
    fractions = phase_fractions(phase)
    assert np.isclose(fractions["crystal"], 0.52)
    assert np.isclose(fractions["oaf"], 0.28)
    assert np.isclose(fractions["maf"], 0.20)
    assert np.isclose(mask.mean(), 0.48)
    assert np.all(mask[phase != CRYSTAL])
    assert not np.any(mask[phase == CRYSTAL])
