from collections import defaultdict
from pathlib import Path

import numpy as np

from pvdf_pf.calibration.bds_full_spectrum import (
    cole_cole_complex,
    fit_relaxation,
    load_figure2_alpha_window,
)
from pvdf_pf.calibration.bds import debye_complex_permittivity
from pvdf_pf.calibration.source_data import load_digitized_curve_csv

ROOT = Path(__file__).resolve().parents[1]
RUI = ROOT / "data/literature/rui2022"


def _static(points, sample, T):
    vals = [p.value for p in points if p.sample_state == sample and p.quantity == "epsilon_c_film" and abs(p.temperature_C-T) < 1e-9]
    assert len(vals) == 1
    return vals[0]


def test_cole_cole_beta_one_reduces_to_debye():
    f = np.logspace(-2, 7, 50)
    cc = cole_cole_complex(f, epsilon_static=20.0, epsilon_infinity=4.0, tau_s=1e-4, beta=1.0)
    debye = debye_complex_permittivity(f, eps_inf=4.0, delta_eps=16.0, tau_s=1e-4)
    assert np.allclose(cc, debye, rtol=1e-12, atol=1e-12)


def test_figure2_digitized_groups_cover_all_same_state_temperatures():
    pts = load_figure2_alpha_window(RUI / "figure_2_bds_alpha_window_digitized.csv")
    counts = defaultdict(int)
    for p in pts:
        counts[(p.sample_state, p.temperature_C)] += 1
        assert p.provenance == "DIGITIZED_SOURCE"
    assert len(counts) == 16
    assert min(counts.values()) >= 10


def test_cole_cole_is_preferred_over_single_debye_for_all_digitized_states():
    pts = load_figure2_alpha_window(RUI / "figure_2_bds_alpha_window_digitized.csv")
    film = load_digitized_curve_csv(RUI / "figure_3a_film_eps_digitized.csv")
    grouped = defaultdict(list)
    for p in pts:
        grouped[(p.sample_state, p.temperature_C)].append(p)
    for (sample, T), group in grouped.items():
        eps_s = _static(film, sample, T)
        deb = fit_relaxation(group, epsilon_static_film=eps_s, model="debye")
        cc = fit_relaxation(group, epsilon_static_film=eps_s, model="cole-cole")
        assert deb.aic - cc.aic > 40.0
        assert 0.20 < cc.beta_cole_cole < 0.60
        assert cc.response_assignment == "combined_OAF_plus_IAF_amorphous_response"
        assert cc.delta_epsilon_amorphous > 0.0


def test_40C_spectra_are_right_censored_by_measurement_ceiling():
    pts = load_figure2_alpha_window(RUI / "figure_2_bds_alpha_window_digitized.csv")
    film = load_digitized_curve_csv(RUI / "figure_3a_film_eps_digitized.csv")
    grouped = defaultdict(list)
    for p in pts:
        grouped[(p.sample_state, p.temperature_C)].append(p)
    for sample in ("unpoled", "poled"):
        cc = fit_relaxation(grouped[(sample, 40.0)], epsilon_static_film=_static(film, sample, 40.0), model="cole-cole")
        assert cc.coverage_status == "RIGHT_CENSORED_PEAK"
        assert cc.peak_coverage_ratio < 3.0
