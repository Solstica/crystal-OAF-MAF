import numpy as np

from pvdf_pf.literature.guo2024 import (
    GUO2024_STRONG_ANISOTROPY_AXIS,
    GUO2024_WEAK_ANISOTROPY_AXIS,
    landau_axis_density,
)
from pvdf_pf.literature.guo2024_vector import (
    GUO2024_STRONG_ANISOTROPY_VECTOR,
    GUO2024_WEAK_ANISOTROPY_VECTOR,
    directional_anisotropy_summary,
    equilibrium_along_direction,
    vector_landau_density,
)


def test_vector_axis_reduces_exactly_to_stage1_polynomial():
    p = np.linspace(-0.12, 0.12, 51)
    for axis_params, vector_params in (
        (GUO2024_STRONG_ANISOTROPY_AXIS, GUO2024_STRONG_ANISOTROPY_VECTOR),
        (GUO2024_WEAK_ANISOTROPY_AXIS, GUO2024_WEAK_ANISOTROPY_VECTOR),
    ):
        vec = np.column_stack((p, np.zeros_like(p), np.zeros_like(p)))
        assert np.allclose(
            vector_landau_density(vec, 25.0, vector_params),
            landau_axis_density(p, 25.0, axis_params),
            rtol=1e-13,
            atol=1e-12,
        )


def test_vector_polynomial_has_sign_and_permutation_symmetry():
    p = np.array([0.011, -0.019, 0.007])
    for params in (
        GUO2024_STRONG_ANISOTROPY_VECTOR,
        GUO2024_WEAK_ANISOTROPY_VECTOR,
    ):
        f0 = vector_landau_density(p, 25.0, params)
        assert np.isclose(f0, vector_landau_density(-p, 25.0, params))
        assert np.isclose(f0, vector_landau_density(p[[2, 0, 1]], 25.0, params))


def test_strong_high_symmetry_directional_minima_regression():
    p100 = equilibrium_along_direction(
        np.array([1.0, 0.0, 0.0]), 25.0, GUO2024_STRONG_ANISOTROPY_VECTOR
    )
    p110 = equilibrium_along_direction(
        np.array([1.0, 1.0, 0.0]), 25.0, GUO2024_STRONG_ANISOTROPY_VECTOR
    )
    p111 = equilibrium_along_direction(
        np.array([1.0, 1.0, 1.0]), 25.0, GUO2024_STRONG_ANISOTROPY_VECTOR
    )
    assert np.isclose(p100["P_abs_min_C_m2"], 0.07301434629379314, rtol=1e-12)
    assert np.isclose(
        p110["P_abs_min_C_m2"],
        np.sqrt(2.0) * p100["P_abs_min_C_m2"],
        rtol=1e-12,
    )
    assert np.isclose(
        p111["P_abs_min_C_m2"],
        np.sqrt(3.0) * p100["P_abs_min_C_m2"],
        rtol=1e-12,
    )
    assert np.isclose(p111["f_min_J_m3"], 3.0 * p100["f_min_J_m3"], rtol=1e-12)


def test_weak_high_symmetry_directional_minima_regression():
    expected = {
        "100": (0.02409590888793086, -949.8339812199513),
        "110": (0.02377743473407198, -953.8016208857639),
        "111": (0.024528720579297194, -1029.4924755961954),
    }
    directions = {
        "100": np.array([1.0, 0.0, 0.0]),
        "110": np.array([1.0, 1.0, 0.0]),
        "111": np.array([1.0, 1.0, 1.0]),
    }
    for key, direction in directions.items():
        result = equilibrium_along_direction(
            direction, 25.0, GUO2024_WEAK_ANISOTROPY_VECTOR
        )
        assert np.isclose(result["P_abs_min_C_m2"], expected[key][0], rtol=1e-12)
        assert np.isclose(result["f_min_J_m3"], expected[key][1], rtol=1e-12)


def test_s16_parameter_sets_separate_strong_from_weak_angular_anisotropy():
    strong = directional_anisotropy_summary(
        25.0, GUO2024_STRONG_ANISOTROPY_VECTOR, n_points=1024
    )
    weak = directional_anisotropy_summary(
        25.0, GUO2024_WEAK_ANISOTROPY_VECTOR, n_points=1024
    )
    assert strong["P_radius_max_over_min"] > 1.70
    assert weak["P_radius_max_over_min"] < 1.04
    assert strong["well_depth_max_over_min"] > 2.9
    assert weak["well_depth_max_over_min"] < 1.09
