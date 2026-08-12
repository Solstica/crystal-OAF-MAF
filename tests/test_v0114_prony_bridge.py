from collections import defaultdict
from pathlib import Path
import math

import numpy as np

from pvdf_pf.calibration.bds_full_spectrum import load_figure2_alpha_window, fit_relaxation
from pvdf_pf.calibration.prony import (
    fit_positive_prony_from_cole_cole,
    generalized_debye_permittivity,
)
from pvdf_pf.calibration.source_data import load_digitized_curve_csv
from pvdf_pf.physics.generalized_debye import (
    GeneralizedDebyeBank,
    advance_generalized_debye_modes,
    continuous_generalized_susceptibility,
    discrete_generalized_susceptibility,
    total_auxiliary_polarization,
)

ROOT = Path(__file__).resolve().parents[1]
RUI = ROOT / "data/literature/rui2022"
TEMPERATURES = (-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0)


def _source_value(points, sample, T):
    vals = [
        p.value for p in points
        if p.sample_state.lower() == sample
        and p.quantity == "epsilon_c_film"
        and abs(p.temperature_C - T) < 1e-9
    ]
    assert len(vals) == 1
    return float(vals[0])


def _all_prony_representations():
    spectra = load_figure2_alpha_window(RUI / "figure_2_bds_alpha_window_digitized.csv")
    static = load_digitized_curve_csv(RUI / "figure_3a_film_eps_digitized.csv")
    grouped = defaultdict(list)
    for point in spectra:
        grouped[(point.sample_state, point.temperature_C)].append(point)
    out = []
    for sample in ("unpoled", "poled"):
        for T in TEMPERATURES:
            cc = fit_relaxation(
                grouped[(sample, T)],
                epsilon_static_film=_source_value(static, sample, T),
                model="cole-cole",
            )
            out.append((cc, fit_positive_prony_from_cole_cole(cc)))
    return out


def test_all_source_states_have_passive_accurate_prony_representation():
    for cc, rep in _all_prony_representations():
        strength = np.asarray(rep.delta_epsilon_modes)
        assert rep.n_modes == 49
        assert np.all(strength >= 0.0)
        assert np.isclose(strength.sum(), cc.delta_epsilon_amorphous, rtol=0, atol=1e-10)
        assert rep.normalized_max_complex_error < 5e-4
        assert rep.normalized_rms_complex_error < 1e-4
        assert rep.active_mode_count >= 30


def test_prony_frequency_response_recovers_cole_cole_static_strength():
    cc, rep = _all_prony_representations()[3]
    eps = generalized_debye_permittivity(
        np.array([1e-8, 1e12]),
        epsilon_infinity=rep.epsilon_infinity_amorphous,
        tau_modes_s=rep.tau_modes_s,
        delta_epsilon_modes=rep.delta_epsilon_modes,
    )
    assert np.isclose(eps[0].real, cc.epsilon_static_amorphous, rtol=2e-3)
    assert np.isclose(eps[1].real, cc.epsilon_infinity_amorphous, rtol=2e-3)


def test_generalized_debye_bank_exact_step_and_mask():
    bank = GeneralizedDebyeBank(
        epsilon_infinity=4.0,
        tau_modes_s=(1e-3, 1e-2),
        delta_epsilon_modes=(2.0, 3.0),
    )
    E = np.ones((2, 3)) * 1e5
    P0 = np.zeros((2, 2, 3))
    mask = np.array([[True, False, True], [False, True, False]])
    P1 = advance_generalized_debye_modes(P0, E, dt_s=2e-3, bank=bank, mask=mask)
    assert P1.shape == P0.shape
    assert np.all(P1[:, ~mask] == 0.0)
    assert np.all(P1[:, mask] > 0.0)
    total = total_auxiliary_polarization(P1)
    assert total.shape == E.shape
    assert np.all(total[mask] > 0.0)


def test_discrete_bank_converges_to_continuous_harmonic_response():
    cc, rep = _all_prony_representations()[8]
    bank = GeneralizedDebyeBank(
        epsilon_infinity=rep.epsilon_infinity_amorphous,
        tau_modes_s=rep.tau_modes_s,
        delta_epsilon_modes=rep.delta_epsilon_modes,
    )
    f = 1.0 / (2.0 * math.pi * rep.target_tau_s)
    continuous = complex(continuous_generalized_susceptibility(bank, f))
    dt = 1.0 / (5000.0 * f)
    discrete = discrete_generalized_susceptibility(bank, frequency_Hz=f, dt_s=dt)
    assert abs(discrete - continuous) / rep.delta_epsilon_total < 5e-4
