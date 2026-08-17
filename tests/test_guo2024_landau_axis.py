import numpy as np

from pvdf_pf.literature.guo2024 import (
    GUO2024_STRONG_ANISOTROPY_AXIS,
    GUO2024_WEAK_ANISOTROPY_AXIS,
    axis_equilibrium,
    landau_axis_density,
    landau_axis_derivative,
)


def test_source_alpha1_law_at_25C():
    expected = -2.4004e6
    assert np.isclose(GUO2024_STRONG_ANISOTROPY_AXIS.alpha1(25.0), expected)
    assert np.isclose(GUO2024_WEAK_ANISOTROPY_AXIS.alpha1(25.0), expected)


def test_axis_polynomial_is_even_and_derivative_is_odd():
    p = np.linspace(0.0, 0.10, 21)
    for params in (GUO2024_STRONG_ANISOTROPY_AXIS, GUO2024_WEAK_ANISOTROPY_AXIS):
        assert np.allclose(
            landau_axis_density(p, 25.0, params),
            landau_axis_density(-p, 25.0, params),
        )
        assert np.allclose(
            landau_axis_derivative(p, 25.0, params),
            -landau_axis_derivative(-p, 25.0, params),
        )


def test_strong_anisotropy_axis_minimum_regression():
    result = axis_equilibrium(25.0, GUO2024_STRONG_ANISOTROPY_AXIS)
    assert np.isclose(result["P_abs_min_C_m2"], 0.07301434629379314, rtol=1e-12)
    assert np.isclose(result["f_min_J_m3"], -25981.404082461617, rtol=1e-12)
    assert abs(
        landau_axis_derivative(
            result["P_abs_min_C_m2"], 25.0, GUO2024_STRONG_ANISOTROPY_AXIS
        )
    ) < 1e-7
    assert result["curvature_at_min_J_m5_C2"] > 0.0


def test_weak_anisotropy_axis_minimum_regression():
    result = axis_equilibrium(25.0, GUO2024_WEAK_ANISOTROPY_AXIS)
    assert np.isclose(result["P_abs_min_C_m2"], 0.02409590888793086, rtol=1e-12)
    assert np.isclose(result["f_min_J_m3"], -949.8339812199513, rtol=1e-12)
    assert abs(
        landau_axis_derivative(
            result["P_abs_min_C_m2"], 25.0, GUO2024_WEAK_ANISOTROPY_AXIS
        )
    ) < 1e-7
    assert result["curvature_at_min_J_m5_C2"] > 0.0


def test_no_cross_coefficient_is_needed_for_axis_slice():
    # This test encodes the scope boundary: only coefficients printed in S2/S3
    # that survive when two polarization components are exactly zero are used.
    for params in (GUO2024_STRONG_ANISOTROPY_AXIS, GUO2024_WEAK_ANISOTROPY_AXIS):
        assert hasattr(params, "alpha11_J_m5_C4")
        assert hasattr(params, "alpha111_J_m9_C6")
        assert not hasattr(params, "alpha12_J_m5_C4")
        assert not hasattr(params, "alpha112_J_m9_C6")
        assert not hasattr(params, "alpha123_J_m9_C6")
