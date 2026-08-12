from pathlib import Path

import numpy as np

from pvdf_pf.calibration.film_closure import closure_series, summarize_closure
from pvdf_pf.calibration.iaf_permittivity import (
    epsilon_iaf_from_melt_extrapolation,
    fit_molten_pvdf_permittivity,
)

ROOT = Path(__file__).resolve().parents[1]
RUI = ROOT / "data/literature/rui2022"
FIG1B = RUI / "figure_1b_melt_eps_digitized.csv"


def test_kirkwood_frohlich_extrapolation_reproduces_rui2021_anchor_values():
    fit = fit_molten_pvdf_permittivity(FIG1B)
    # Rui 2021 reports epsilon_am ~= 19.5 at 25 C and 12.1 at 127 C.
    assert np.isclose(fit.epsilon(25.0), 19.5, atol=0.15)
    assert np.isclose(fit.epsilon(127.0), 12.1, atol=0.15)
    assert fit.g_r_squared > 0.99


def test_kirkwood_frohlich_iaf_is_not_direct_linear_epsilon_extrapolation():
    cold, _ = epsilon_iaf_from_melt_extrapolation(FIG1B, temperature_C=-30.0)
    roomish, _ = epsilon_iaf_from_melt_extrapolation(FIG1B, temperature_C=20.0)
    hot, _ = epsilon_iaf_from_melt_extrapolation(FIG1B, temperature_C=40.0)
    assert 25.5 < cold < 26.8
    assert 19.5 < roomish < 20.6
    assert 17.7 < hot < 18.6
    assert cold > roomish > hot


def test_published_s2_values_are_executable_with_nonclosing_raw_mobility_weights():
    for sample, rmse_limit in (("unpoled", 0.10), ("poled", 0.20)):
        rows = closure_series(RUI, sample_state=sample)
        summary = summarize_closure(rows)
        assert summary["rmse_source_reconstruction"] < rmse_limit
        assert summary["max_abs_relative_error_source_reconstruction"] < 0.025
        # This is the essential guardrail: the numerical source reconstruction
        # is not a valid phase-fraction map because the weights sum to ~1.21.
        assert 1.205 < summary["mean_source_reconstruction_fraction_sum"] < 1.215
        assert all(1.205 < r.source_reconstruction_fraction_sum < 1.215 for r in rows)
