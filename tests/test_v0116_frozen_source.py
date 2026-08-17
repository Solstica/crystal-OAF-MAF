import numpy as np

from pvdf_pf.core.grid import Grid2D
from pvdf_pf.morphology.oriented import periodic_winding_three_phase
from pvdf_pf.morphology.three_phase import CRYSTAL, OAF
from pvdf_pf.physics.coupled_local_field import (
    advance_self_consistent_generalized_debye,
    build_combined_amorphous_background_map,
)
from pvdf_pf.physics.frozen_source import (
    advance_generalized_debye_with_frozen_source,
    build_lamellar_frozen_polarization_z,
)
from pvdf_pf.physics.generalized_debye import GeneralizedDebyeBank


HUANG_FC = 5.78 / 11.8
HUANG_FO = 3.02 / 11.8
HUANG_FI = 3.00 / 11.8


def _bank():
    return GeneralizedDebyeBank(
        epsilon_infinity=5.0,
        tau_modes_s=(0.2, 1.5),
        delta_epsilon_modes=(1.0, 2.0),
    )


def test_zero_frozen_source_reproduces_v0115_step():
    grid = Grid2D(nz=24, nx=20)
    phase = periodic_winding_three_phase(
        grid,
        crystal_fraction=0.50,
        oaf_fraction=0.25,
        maf_fraction=0.25,
        winding_z=1,
        winding_x=0,
    )
    bank = _bank()
    eps_b, mask = build_combined_amorphous_background_map(
        phase, bank, epsilon_crystal_background=3.0
    )
    P0 = np.zeros((bank.n_modes,) + grid.shape)
    kwargs = dict(E_external_z=8.0e5, dt_s=0.3, bank=bank)

    baseline = advance_self_consistent_generalized_debye(
        P0, eps_b, mask, grid, **kwargs
    )
    frozen = advance_generalized_debye_with_frozen_source(
        P0, eps_b, mask, np.zeros(grid.shape), grid, **kwargs
    )

    assert np.allclose(frozen.E_z, baseline.E_z, rtol=1e-12, atol=1e-10)
    assert np.allclose(frozen.P_modes, baseline.P_modes, rtol=1e-12, atol=1e-20)
    assert np.isclose(frozen.D_z_mean_face, baseline.D_z_mean_face, rtol=1e-12)
    assert np.isclose(frozen.D_z_mean_cell, baseline.D_z_mean_cell, rtol=1e-12)


def test_huang_aligned_lamellae_tangent_polarization_has_no_bound_charge_field():
    grid = Grid2D(nz=36, nx=36)
    phase = periodic_winding_three_phase(
        grid,
        crystal_fraction=HUANG_FC,
        oaf_fraction=HUANG_FO,
        maf_fraction=HUANG_FI,
        winding_z=0,
        winding_x=1,
    )
    bank = _bank()
    eps_b, mask = build_combined_amorphous_background_map(
        phase, bank, epsilon_crystal_background=3.0
    )
    P_frozen = build_lamellar_frozen_polarization_z(
        phase,
        grid,
        crystal_fraction=HUANG_FC,
        oaf_fraction=HUANG_FO,
        p_beta_C_m2=0.010,
        p_oaf_mean_C_m2=0.006,
        winding_z=0,
        winding_x=1,
    )
    P0 = np.zeros((bank.n_modes,) + grid.shape)
    step = advance_generalized_debye_with_frozen_source(
        P0,
        eps_b,
        mask,
        P_frozen,
        grid,
        E_external_z=0.0,
        dt_s=0.4,
        bank=bank,
    )

    assert np.max(np.abs(step.E_z)) < 1e-8
    assert np.max(np.abs(step.E_x)) < 1e-8
    assert np.max(np.abs(step.P_total)) < 1e-18
    assert step.gauss_relative_residual < 1e-10


def test_same_frozen_source_normal_to_layers_generates_internal_field():
    grid = Grid2D(nz=36, nx=36)
    phase = periodic_winding_three_phase(
        grid,
        crystal_fraction=HUANG_FC,
        oaf_fraction=HUANG_FO,
        maf_fraction=HUANG_FI,
        winding_z=1,
        winding_x=0,
    )
    bank = _bank()
    eps_b, mask = build_combined_amorphous_background_map(
        phase, bank, epsilon_crystal_background=3.0
    )
    P_frozen = build_lamellar_frozen_polarization_z(
        phase,
        grid,
        crystal_fraction=HUANG_FC,
        oaf_fraction=HUANG_FO,
        p_beta_C_m2=0.010,
        p_oaf_mean_C_m2=0.006,
        winding_z=1,
        winding_x=0,
    )
    P0 = np.zeros((bank.n_modes,) + grid.shape)
    step = advance_generalized_debye_with_frozen_source(
        P0,
        eps_b,
        mask,
        P_frozen,
        grid,
        E_external_z=0.0,
        dt_s=0.4,
        bank=bank,
    )

    assert np.max(np.abs(step.E_z)) > 1.0e6
    assert abs(float(np.mean(step.E_z))) < 1e-7 * np.max(np.abs(step.E_z))
    assert step.gauss_relative_residual < 1e-7
    # The solver enforces face-flux Gauss balance.  The diagnostic below compares
    # that face flux with a cell-centered constitutive reconstruction, so a small
    # O(dx) mismatch remains at discontinuous beta/OAF/IAF interfaces.
    assert step.constitutive_flux_mismatch < 1e-3


def test_interface_decay_preserves_discrete_oaf_mean():
    grid = Grid2D(nz=48, nx=48)
    phase = periodic_winding_three_phase(
        grid,
        crystal_fraction=HUANG_FC,
        oaf_fraction=HUANG_FO,
        maf_fraction=HUANG_FI,
        winding_z=1,
        winding_x=0,
    )
    target = 0.007
    P_uniform = build_lamellar_frozen_polarization_z(
        phase,
        grid,
        crystal_fraction=HUANG_FC,
        oaf_fraction=HUANG_FO,
        p_beta_C_m2=0.010,
        p_oaf_mean_C_m2=target,
        winding_z=1,
        winding_x=0,
        oaf_profile="uniform",
    )
    P_decay = build_lamellar_frozen_polarization_z(
        phase,
        grid,
        crystal_fraction=HUANG_FC,
        oaf_fraction=HUANG_FO,
        p_beta_C_m2=0.010,
        p_oaf_mean_C_m2=target,
        winding_z=1,
        winding_x=0,
        oaf_profile="interface_decay",
        decay_length_fraction_of_half_oaf=0.35,
    )
    oaf = phase == OAF
    crystal = phase == CRYSTAL

    assert np.isclose(np.mean(P_uniform[oaf]), target, rtol=0, atol=1e-15)
    assert np.isclose(np.mean(P_decay[oaf]), target, rtol=0, atol=1e-15)
    assert np.max(P_decay[oaf]) > target
    assert np.min(P_decay[oaf]) < target
    assert np.allclose(P_uniform[crystal], 0.010)
    assert np.allclose(P_decay[crystal], 0.010)
