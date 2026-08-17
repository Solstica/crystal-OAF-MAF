import numpy as np

from pvdf_pf.literature.ahluwalia2008_dimensionless import (
    dimensionless_checkpoint_300K,
    inferred_length_from_noise_z,
    polarization_scale,
)


def test_source_polarization_scale_regression():
    assert np.isclose(polarization_scale(), 0.10178585540417863, rtol=1e-13)


def test_reported_300K_noise_implies_subnanometer_rescaling_length():
    xi = inferred_length_from_noise_z(300.0)
    assert np.isclose(xi, 4.0481614607590865e-10, rtol=1e-12)
    assert np.isclose(xi * 1e9, 0.40481614607590865, rtol=1e-12)


def test_inferred_rescaling_length_is_not_the_2p16nm_grid_spacing():
    c = dimensionless_checkpoint_300K()
    assert c.physical_grid_spacing_nm == 2.16
    assert c.inferred_length_scale_nm < 0.5
    assert c.grid_spacing_dimensionless > 5.0
    assert np.isclose(c.grid_spacing_dimensionless, 5.335755554559748, rtol=1e-12)


def test_300K_dimensionless_coefficients_regression():
    c = dimensionless_checkpoint_300K()
    assert np.isclose(c.alpha_prime, 0.04542396975425331, rtol=1e-12)
    assert np.isclose(c.transverse_susceptibility_m_F, 6.089418029639882e-12, rtol=1e-12)
    assert np.isclose(c.alpha_xx_prime, 1.7229057576070796, rtol=1e-12)
    assert np.isclose(c.alpha_yy_prime, c.alpha_xx_prime, rtol=1e-14)
    assert np.isclose(c.E_prime, 1.1849187757883282, rtol=1e-12)
    assert np.isclose(c.K1_prime, 1.3495602114497969, rtol=1e-12)
    assert np.isclose(c.K2_prime, c.K1_prime, rtol=1e-14)
    assert np.isclose(c.K3_prime, c.K1_prime, rtol=1e-14)


def test_300K_state_and_kinetic_ratios_regression():
    c = dimensionless_checkpoint_300K()
    assert np.isclose(c.equilibrium_Pr_C_m2, 0.09932864465590406, rtol=1e-12)
    assert np.isclose(c.equilibrium_w_dimensionless, 0.9758590155919278, rtol=1e-12)
    assert np.isclose(c.initial_w_dimensionless, 1.1199984472038926, rtol=1e-12)
    assert np.isclose(c.m_Gamma_x_over_Gamma_z, 0.0399775059754201, rtol=1e-12)
    assert np.isclose(c.n_Gamma_y_over_Gamma_z, 9.001687368170652, rtol=1e-12)
    assert c.ps_per_tstar == 9.0
    assert np.isclose(c.cell_length_dimensionless, 341.4883554918239, rtol=1e-12)
