import numpy as np

from pvdf_pf.calibration.phase_fraction_state import reconstruct_above_tg
from pvdf_pf.morphology.oaf_devitrification import (
    SOURCE_CRYSTAL,
    SOURCE_IAF,
    SOURCE_MOAF,
    SOURCE_ROAF,
    allocate_source_oaf_subphases,
    source_phase_fractions,
)
from pvdf_pf.morphology.three_phase import CRYSTAL, MAF, OAF


def _base_map():
    # 40 pixels: 24 crystal, 8 OAF, 8 project-MAF -> source 0.6/0.2/0.2.
    phase = np.array(
        [CRYSTAL] * 24 + [OAF] * 8 + [MAF] * 8,
        dtype=np.int8,
    ).reshape(5, 8)
    return phase


def test_spatial_allocation_hits_fraction_to_pixel_resolution():
    phase = _base_map()
    state = reconstruct_above_tg(
        x_raf=0.10,
        x_maf=0.30,
        temperature_C=20.0,
        sample_state="poled",
    )
    # MOAF = 0.10, total OAF = 0.20 -> 50% of OAF pixels mobile.
    score = np.arange(phase.size, dtype=float).reshape(phase.shape)
    source, report = allocate_source_oaf_subphases(
        phase, fraction_state=state, mobility_score=score
    )
    assert report.n_oaf_pixels == 8
    assert report.n_moaf_pixels == 4
    assert np.isclose(report.realized_mobile_fraction_within_oaf, 0.5)
    assert np.count_nonzero(source == SOURCE_CRYSTAL) == 24
    assert np.count_nonzero(source == SOURCE_ROAF) == 4
    assert np.count_nonzero(source == SOURCE_MOAF) == 4
    assert np.count_nonzero(source == SOURCE_IAF) == 8
    fractions = source_phase_fractions(source)
    assert np.isclose(fractions["crystal"], 0.6)
    assert np.isclose(fractions["roaf"], 0.1)
    assert np.isclose(fractions["moaf"], 0.1)
    assert np.isclose(fractions["iaf"], 0.2)


def test_highest_scores_devitrify_first():
    phase = _base_map()
    state = reconstruct_above_tg(
        x_raf=0.15,
        x_maf=0.25,
        temperature_C=10.0,
        sample_state="unpoled",
    )
    # MOAF = 0.05 => 25% of OAF = two of eight OAF pixels.
    score = np.zeros(phase.shape, dtype=float)
    oaf_idx = np.flatnonzero(phase.ravel() == OAF)
    score.ravel()[oaf_idx] = np.arange(1, 9, dtype=float)
    source, report = allocate_source_oaf_subphases(
        phase, fraction_state=state, mobility_score=score
    )
    selected = np.flatnonzero(source.ravel() == SOURCE_MOAF)
    assert set(selected) == set(oaf_idx[-2:])
    assert report.score_min_selected == 7.0
    assert report.score_max_unselected == 6.0
