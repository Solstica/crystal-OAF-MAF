import numpy as np

from pvdf_pf.calibration.bds import (
    debye_complex_permittivity,
    effective_dipole_moment_from_g,
    fit_debye_spectrum,
    kirkwood_frohlich_g_mu2,
    kirkwood_g_from_mu,
)


def test_debye_fit_recovers_synthetic_parameters():
    frequency = np.logspace(0, 5, 121)
    expected = {"eps_inf": 4.2, "delta_eps": 17.5, "tau_s": 2.5e-3}
    eps = debye_complex_permittivity(frequency, **expected)
    fit = fit_debye_spectrum(frequency, eps.real, -eps.imag)

    assert fit.success
    assert np.isclose(fit.eps_inf, expected["eps_inf"], rtol=1e-5)
    assert np.isclose(fit.delta_eps, expected["delta_eps"], rtol=1e-5)
    assert np.isclose(fit.tau_s, expected["tau_s"], rtol=1e-5)
    assert fit.rmse_real < 1e-6
    assert fit.rmse_loss < 1e-6


def test_kirkwood_inversions_are_algebraically_consistent():
    kwargs = {
        "epsilon_static": 18.0,
        "epsilon_fast": 4.0,
        "active_dipole_number_density_m3": 1.2e28,
        "temperature_K": 300.0,
    }
    g_mu2 = kirkwood_frohlich_g_mu2(**kwargs)
    mu = 7.0e-30
    g = kirkwood_g_from_mu(**kwargs, dipole_moment_Cm=mu)
    mu_back = effective_dipole_moment_from_g(**kwargs, kirkwood_g=g)

    assert g_mu2 > 0.0
    assert g > 0.0
    assert np.isclose(mu_back, mu, rtol=1e-12)
