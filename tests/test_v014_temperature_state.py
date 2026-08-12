from pathlib import Path

import numpy as np
import pytest

from pvdf_pf.calibration.rui2022_state import (
    build_rui2022_temperature_state,
    interpolate_digitized_curve,
)
from pvdf_pf.calibration.source_data import load_digitized_curve_csv


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "literature" / "rui2022"


def test_exact_source_point_is_reproduced():
    state = build_rui2022_temperature_state(
        DATA, temperature_C=0.0, sample_state="unpoled", mobility_boundary="l1"
    )
    assert np.isclose(state.epsilon_moaf, 28.66, atol=1e-12)
    assert np.isclose(state.interaction_lambda, 1.5, atol=1e-12)
    assert np.isclose(state.rotational_mobility_m2_per_Vs, 7.07e5, rtol=1e-12)


def test_linear_interpolation_for_n():
    s1 = load_digitized_curve_csv(DATA / "figure_s1_digitized.csv")
    value = interpolate_digitized_curve(
        s1,
        -25.0,
        quantity="n",
        sample_state="unpoled",
        panel="A",
    )
    assert np.isclose(value, 0.5 * (2.695 + 3.08), atol=1e-12)


def test_mu_r_uses_log_interpolation():
    state = build_rui2022_temperature_state(
        DATA, temperature_C=-25.0, sample_state="unpoled", mobility_boundary="l1"
    )
    expected = np.sqrt(1.11e3 * 1.69e4)
    assert np.isclose(state.rotational_mobility_m2_per_Vs, expected, rtol=1e-12)


def test_no_extrapolation():
    with pytest.raises(ValueError, match="outside source range"):
        build_rui2022_temperature_state(
            DATA, temperature_C=45.0, sample_state="poled", mobility_boundary="l2"
        )


def test_devitrification_index_is_diagnostic_endpoints():
    low = build_rui2022_temperature_state(
        DATA, temperature_C=-30.0, sample_state="poled", mobility_boundary="l1"
    )
    high = build_rui2022_temperature_state(
        DATA, temperature_C=40.0, sample_state="poled", mobility_boundary="l1"
    )
    assert np.isclose(low.active_dipole_devitrification_index, 0.0)
    assert np.isclose(high.active_dipole_devitrification_index, 1.0)
    assert high.rotational_mobility_decades_from_minus30C > 4.0
