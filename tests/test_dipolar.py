import numpy as np

from pvdf_pf.physics.dipolar import (
    EPS0,
    advance_debye_polarization,
    complex_susceptibility_from_harmonic_history,
    simulate_debye_history,
)


def test_exact_step_response():
    delta_eps = 12.0
    tau = 2.0e-3
    dt = 1.0e-3
    E = 2.0e6
    target = EPS0 * delta_eps * E

    P1 = advance_debye_polarization(0.0, E, dt_s=dt, tau_s=tau, delta_eps=delta_eps)
    expected = target * (1.0 - np.exp(-dt / tau))
    assert np.isclose(P1, expected, rtol=1e-14)


def test_oaf_mask_zeroes_non_oaf_cells():
    P = np.zeros((2, 2))
    E = np.ones((2, 2)) * 1e6
    mask = np.array([[True, False], [False, True]])
    out = advance_debye_polarization(
        P,
        E,
        dt_s=1e-4,
        tau_s=1e-3,
        delta_eps=10.0,
        mask=mask,
    )
    assert np.all(out[~mask] == 0.0)
    assert np.all(out[mask] > 0.0)


def test_harmonic_time_domain_matches_debye_frequency_response():
    frequency = 25.0
    tau = 3.0e-3
    delta_eps = 9.0
    steps_per_cycle = 400
    cycles = 40
    dt = 1.0 / (frequency * steps_per_cycle)
    t = dt * np.arange(steps_per_cycle * cycles)
    E = 1.5e6 * np.cos(2.0 * np.pi * frequency * t)
    P = simulate_debye_history(E, dt_s=dt, tau_s=tau, delta_eps=delta_eps)
    chi = complex_susceptibility_from_harmonic_history(
        E,
        P,
        dt_s=dt,
        frequency_hz=frequency,
        discard_fraction=0.5,
    )
    expected = delta_eps / (1.0 + 1j * 2.0 * np.pi * frequency * tau)
    assert np.isclose(chi.real, expected.real, rtol=3e-3, atol=3e-3)
    assert np.isclose(chi.imag, expected.imag, rtol=3e-3, atol=3e-3)
