from pathlib import Path

import numpy as np

from pvdf_pf.calibration.source_state_transfer import (
    closed_oaf_series,
    source_state_definitions,
    summarize_closed_oaf,
)

ROOT = Path(__file__).resolve().parents[1]
RUI = ROOT / "data/literature/rui2022"


def test_source_states_make_sample_transfer_explicit():
    states = {s.name: s for s in source_state_definitions()}
    assert np.isclose(states["rui2022_bds_target"].crystal_fraction, 0.52)
    assert np.isclose(states["rui2022_tmdsc_mobility_donor"].crystal_fraction, 0.59)
    assert np.isclose(states["rui2022_si_s2_rounded_model"].crystal_fraction, 0.60)
    assert np.isclose(states["rui2021_poled_bopvdf_structural_proxy"].closed_sum, 1.0)
    assert states["rui2022_bds_target"].sample_state != states["rui2022_tmdsc_mobility_donor"].sample_state


def test_same_state_poled_oaf_effective_response_is_positive_and_increases_over_window():
    rows = closed_oaf_series(RUI, sample_state="poled")
    eps = np.asarray([r.epsilon_oaf_effective_required for r in rows])
    assert np.all(eps > 0.0)
    assert eps[-1] > eps[0]
    assert 25.0 < eps[0] < 30.0
    assert 44.0 < eps[-1] < 48.0


def test_same_state_unpoled_proxy_is_reported_only_as_comparison_but_numerically_sensible():
    rows = closed_oaf_series(RUI, sample_state="unpoled")
    eps = np.asarray([r.epsilon_oaf_effective_required for r in rows])
    assert np.all(eps > 0.0)
    assert 16.0 < eps[0] < 21.0
    assert 28.0 < eps[-1] < 33.0


def test_effective_oaf_and_figure_s2_moaf_are_not_conflated():
    rows = closed_oaf_series(RUI, sample_state="poled")
    summary = summarize_closed_oaf(rows)
    assert summary["rmse_vs_figure_s2_moaf"] > 5.0
    assert rows[0].epsilon_oaf_effective_required < rows[0].epsilon_moaf_figure_s2
    assert rows[-1].epsilon_oaf_effective_required > rows[-1].epsilon_moaf_figure_s2
