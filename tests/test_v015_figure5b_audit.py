from pathlib import Path

import numpy as np

from pvdf_pf.calibration.figure5b_audit import audit_figure5b, source_anchor_checks


ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "data/literature/rui2022/figure_5b_maintext_digitized.csv"


def test_main_text_anchors_are_preserved():
    checks = source_anchor_checks(CSV)
    assert checks
    assert all(checks.values())


def test_figure5b_closes_to_raw_amorphous_fraction():
    points = audit_figure5b(CSV)
    residuals = np.asarray([p.raw_closure_residual_vs_0p41 for p in points])
    assert np.max(np.abs(residuals)) <= 0.012


def test_literal_si_mapping_exposes_near_tg_inconsistency():
    points = audit_figure5b(CSV)
    p_minus30 = next(p for p in points if np.isclose(p.temperature_C, -30.0))
    assert np.isclose(p_minus30.x_maf, 0.166)
    assert np.isclose(p_minus30.si_implied_eta_moaf, -0.034)
    assert not p_minus30.si_nonnegative
    assert not p_minus30.literal_si_mapping_valid


def test_high_temperature_point_is_nonnegative_but_retains_rounding_residual():
    points = audit_figure5b(CSV)
    p40 = next(p for p in points if np.isclose(p.temperature_C, 40.0))
    assert np.isclose(p40.si_implied_eta_moaf, 0.196)
    assert p40.si_nonnegative
    assert np.isclose(p40.si_closure_residual, 0.010, atol=1e-6)
