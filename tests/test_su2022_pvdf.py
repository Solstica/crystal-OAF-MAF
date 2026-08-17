import numpy as np

from pvdf_pf.literature.su2022 import (
    SU2022_PVDF,
    electrostatic_displacement,
    pvdf_bulk_density,
    z_axis_global_minimum,
    z_axis_stationary_points,
)


def test_source_table_3_coefficients_are_transcribed():
    assert SU2022_PVDF.alpha1_J_m_C2 == 5.647e9
    assert SU2022_PVDF.alpha2_J_m_C2 == 5.647e9
    assert SU2022_PVDF.alpha3_prefactor_J_m_C2_K == 1.412e7
    assert SU2022_PVDF.alpha3_zero_K == 315.0
    assert SU2022_PVDF.alpha33_J_m5_C4 == -1.842e11
    assert SU2022_PVDF.alpha333_J_m9_C6 == 2.585e13
    assert SU2022_PVDF.Q11_m4_C2 == -8.5
    assert SU2022_PVDF.Q12_m4_C2 == 0.0
    assert SU2022_PVDF.Q44_m4_C2 == 0.0
    assert SU2022_PVDF.s11_m2_N == 4.0e-10
    assert SU2022_PVDF.s12_m2_N == 1.11e-9
    assert SU2022_PVDF.s44_m2_N == 1.25e-9


def test_alpha3_zero_is_exactly_315_K():
    assert SU2022_PVDF.alpha3(315.0) == 0.0
    assert SU2022_PVDF.alpha3(300.0) == -2.118e8


def test_eq4_has_paraelectric_transverse_and_uniaxial_higher_order_terms():
    t = 300.0
    p = 0.04
    assert np.isclose(
        pvdf_bulk_density(np.array([p, 0.0, 0.0]), t),
        SU2022_PVDF.alpha1_J_m_C2 * p**2,
    )
    expected_z = (
        SU2022_PVDF.alpha3(t) * p**2
        + SU2022_PVDF.alpha33_J_m5_C4 * p**4
        + SU2022_PVDF.alpha333_J_m9_C6 * p**6
    )
    assert np.isclose(pvdf_bulk_density(np.array([0.0, 0.0, p]), t), expected_z)


def test_bulk_energy_is_even_in_all_polarization_components():
    p = np.array([0.011, -0.021, 0.057])
    f = pvdf_bulk_density(p, 300.0)
    assert np.isclose(f, pvdf_bulk_density(-p, 300.0))
    assert np.isclose(f, pvdf_bulk_density(p * np.array([-1.0, 1.0, -1.0]), 300.0))


def test_declared_300K_algebraic_checkpoint_regression():
    result = z_axis_global_minimum(300.0)
    assert result["kind"] == "minimum"
    assert np.isclose(result["Pz_C_m2"], 0.0725867804739834, rtol=1e-12)
    assert np.isclose(result["f_J_m3"], -2448466.201076933, rtol=1e-12)

    points = z_axis_stationary_points(300.0)
    assert points[0]["Pz_C_m2"] == 0.0
    assert points[0]["kind"] == "maximum"


def test_electrostatic_equation_keeps_epsb_explicit():
    e = np.array([1.0e5, 0.0, -2.0e5])
    p = np.array([0.01, 0.0, 0.02])
    d = electrostatic_displacement(e, p, eps_b=10.0)
    expected = 8.8541878128e-12 * 10.0 * e + p
    assert np.allclose(d, expected)
