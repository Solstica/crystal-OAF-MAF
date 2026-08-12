from pathlib import Path

import numpy as np

from pvdf_pf.calibration.same_state_bds import (
    amorphous_complex_permittivity_debye_bridge,
    build_same_state_bds_point,
    film_complex_permittivity_parallel_bridge,
    fit_high_temperature_arrhenius,
    same_state_bds_series,
)

ROOT = Path(__file__).resolve().parents[1]
RUI = ROOT / "data/literature/rui2022"


def test_same_state_two_phase_inversion_closes_film_static_permittivity():
    point = build_same_state_bds_point(
        RUI,
        sample_state="poled",
        temperature_C=20.0,
    )
    reconstructed = (
        point.eta_crystal * point.epsilon_crystal
        + (1.0 - point.eta_crystal) * point.epsilon_amorphous_static
    )
    assert np.isclose(reconstructed, point.epsilon_film_static, atol=1e-12)
    assert point.response_assignment == "combined_OAF_plus_IAF_amorphous_response"


def test_poled_combined_amorphous_static_response_exceeds_unpoled_at_all_points():
    unpoled = same_state_bds_series(RUI, sample_state="unpoled")
    poled = same_state_bds_series(RUI, sample_state="poled")
    assert len(unpoled) == len(poled) == 8
    assert all(
        p.epsilon_amorphous_static > u.epsilon_amorphous_static
        for u, p in zip(unpoled, poled)
    )


def test_digitized_figure3b_recovers_reported_arrhenius_scale():
    reported = {"unpoled": 57.8, "poled": 68.3}
    for sample, target in reported.items():
        fit = fit_high_temperature_arrhenius(
            same_state_bds_series(RUI, sample_state=sample)
        )
        # Figure digitization is approximate; this test checks scale/point identity,
        # not equality to the article's stronger text-reported fit.
        assert fit.r_squared > 0.95
        assert abs(fit.activation_energy_kj_mol - target) / target < 0.15


def test_alpha_relaxation_time_decreases_with_temperature_for_each_film_state():
    for sample in ("unpoled", "poled"):
        rows = same_state_bds_series(RUI, sample_state=sample)
        tau = np.asarray([p.tau_alpha_s for p in rows])
        assert np.all(np.diff(tau) < 0.0)


def test_debye_bridge_requires_external_fast_limit_and_has_correct_limits():
    point = build_same_state_bds_point(
        RUI,
        sample_state="poled",
        temperature_C=20.0,
    )
    eps_fast = 8.0
    eps_am = amorphous_complex_permittivity_debye_bridge(
        np.asarray([0.0, 1e30]),
        point,
        epsilon_amorphous_fast=eps_fast,
    )
    assert np.isclose(eps_am[0].real, point.epsilon_amorphous_static, rtol=1e-12)
    assert np.isclose(eps_am[-1].real, eps_fast, rtol=1e-12)

    eps_film = film_complex_permittivity_parallel_bridge(
        np.asarray([0.0]),
        point,
        epsilon_amorphous_fast=eps_fast,
    )
    assert np.isclose(eps_film[0].real, point.epsilon_film_static, rtol=1e-12)


def test_v0112_does_not_infer_oaf_only_response():
    point = build_same_state_bds_point(
        RUI,
        sample_state="unpoled",
        temperature_C=0.0,
    )
    assert "OAF_plus_IAF" in point.response_assignment
    assert not hasattr(point, "epsilon_oaf")
    assert not hasattr(point, "tau_oaf_s")
