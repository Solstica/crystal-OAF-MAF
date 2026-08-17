import numpy as np

from pvdf_pf.literature.ahluwalia2008 import (
    GRADIENT_CHECKPOINT,
    KINETIC_CHECKPOINT,
    TABLE_I_PRINTED,
    TABLE_II_PRINTED,
    derive_lgd_from_table_i,
    kinetic_ratios,
    table_ii_rounding_audit,
)


def test_printed_source_tables_are_transcribed():
    assert TABLE_I_PRINTED.P0_C_m2 == 0.111
    assert TABLE_I_PRINTED.Tc_K == 450.0
    assert TABLE_I_PRINTED.Pc_C_m2 == 0.088
    assert TABLE_I_PRINTED.f0_J_m3 == -2.12e8

    assert TABLE_II_PRINTED.alpha0_J_m_C2_K == 9.02e7
    assert TABLE_II_PRINTED.T0_K == 252.0
    assert TABLE_II_PRINTED.beta_J_m5_C4 == 9.20e12
    assert TABLE_II_PRINTED.gamma_J_m9_C6 == 8.88e14


def test_eqs_2_3_applied_to_rounded_table_i_are_close_but_not_identical_to_table_ii():
    derived = derive_lgd_from_table_i()
    assert np.isclose(derived.alpha0_J_m_C2_K, 9.067692423802367e7, rtol=1e-12)
    assert np.isclose(derived.T0_K, 248.21125396728488, rtol=1e-12)
    assert np.isclose(derived.beta_J_m5_C4, 9.451230803767729e12, rtol=1e-12)
    assert np.isclose(derived.gamma_J_m9_C6, 9.153438924103561e14, rtol=1e-12)

    audit = table_ii_rounding_audit()
    differences = [abs(item["relative_difference"]) for item in audit.values()]
    assert max(differences) < 0.04
    assert max(differences) > 0.02


def test_md_domain_wall_checkpoint_and_k3_scope_are_explicit():
    assert GRADIENT_CHECKPOINT.xi1_m == 0.4e-9
    assert GRADIENT_CHECKPOINT.xi2_m == 0.4e-9
    assert GRADIENT_CHECKPOINT.K1_J_m3_C2 == 2.108e-8
    assert GRADIENT_CHECKPOINT.K2_J_m3_C2 == 2.108e-8
    assert GRADIENT_CHECKPOINT.K3_J_m3_C2 == 2.108e-8
    assert "CONVENIENCE" in GRADIENT_CHECKPOINT.K3_status


def test_noise_matching_recovers_source_kinetic_ratios_and_time_mapping():
    ratios = kinetic_ratios()
    assert np.isclose(ratios["m_Gamma_x_over_Gamma_z"], 0.04, rtol=1e-3)
    assert np.isclose(ratios["n_Gamma_y_over_Gamma_z"], 9.0, rtol=2e-4)
    assert ratios["ps_per_tstar_from_equilibration_match"] == 9.0

    assert KINETIC_CHECKPOINT.reported_smallest_tdgl_time_ps == 9.0
    assert KINETIC_CHECKPOINT.tdgl_grid_length_nm == 2.16
