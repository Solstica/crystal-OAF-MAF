import numpy as np

from pvdf_pf.core.grid import Grid2D
from pvdf_pf.morphology.oriented import wavy_winding_three_phase
from pvdf_pf.physics.coupled_local_field import build_combined_amorphous_background_map
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


def _solve(amplitude: float):
    grid = Grid2D(nz=48, nx=48)
    phase = wavy_winding_three_phase(
        grid,
        crystal_fraction=HUANG_FC,
        oaf_fraction=HUANG_FO,
        maf_fraction=HUANG_FI,
        winding_z=0,
        winding_x=1,
        waviness_amplitude=amplitude,
        waviness_mode_z=1,
        waviness_mode_x=0,
    )
    bank = _bank()
    eps_b, amorphous = build_combined_amorphous_background_map(
        phase,
        bank,
        epsilon_crystal_background=3.0,
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
        oaf_profile="uniform",
    )
    P0 = np.zeros((bank.n_modes,) + grid.shape)
    step = advance_generalized_debye_with_frozen_source(
        P0,
        eps_b,
        amorphous,
        P_frozen,
        grid,
        E_external_z=0.0,
        dt_s=0.4,
        bank=bank,
    )
    source_z_variation = float(
        np.sqrt(np.mean((P_frozen - np.roll(P_frozen, 1, axis=0)) ** 2))
    )
    return step, source_z_variation


def test_zero_waviness_is_tangent_interface_limit():
    step, source_z_variation = _solve(0.0)
    assert source_z_variation == 0.0
    assert np.max(np.abs(step.E_z)) < 1e-8
    assert np.max(np.abs(step.E_x)) < 1e-8
    assert step.gauss_relative_residual < 1e-10


def test_waviness_creates_bound_charge_field_without_rotating_mean_lamella_orientation():
    step, source_z_variation = _solve(0.10)
    rms_field = float(np.sqrt(np.mean(step.E_z**2 + step.E_x**2)))
    assert source_z_variation > 0.0
    assert rms_field > 1.0e5
    assert np.max(np.abs(step.E_z)) > 1.0e5
    assert step.gauss_relative_residual < 1e-7


def test_larger_waviness_increases_geometry_induced_field_over_small_amplitude_case():
    small, _ = _solve(0.02)
    large, _ = _solve(0.10)
    rms_small = float(np.sqrt(np.mean(small.E_z**2 + small.E_x**2)))
    rms_large = float(np.sqrt(np.mean(large.E_z**2 + large.E_x**2)))
    assert rms_small > 0.0
    assert rms_large > rms_small
