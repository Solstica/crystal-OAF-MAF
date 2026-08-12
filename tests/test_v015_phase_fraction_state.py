from pathlib import Path

import numpy as np
import pytest

from pvdf_pf.calibration.phase_fraction_state import (
    below_tg_limit,
    reconstruct_above_tg,
    reconstruct_curve_from_figure5b,
)


def test_above_tg_reconstruction_closes_source_model():
    state = reconstruct_above_tg(
        x_raf=0.12,
        x_maf=0.28,
        temperature_C=20.0,
        sample_state="poled",
        crystal_fraction=0.6,
        iaf_fraction=0.2,
    )
    assert np.isclose(state.crystal, 0.6)
    assert np.isclose(state.roaf, 0.12)
    assert np.isclose(state.moaf, 0.08)
    assert np.isclose(state.iaf, 0.2)
    assert np.isclose(state.total_oaf, 0.2)
    assert np.isclose(state.mobile_fraction_within_oaf, 0.4)


def test_above_tg_rejects_maf_smaller_than_iaf():
    with pytest.raises(ValueError):
        reconstruct_above_tg(
            x_raf=0.25,
            x_maf=0.15,
            temperature_C=0.0,
            sample_state="unpoled",
        )


def test_below_tg_limit_keeps_oaf_rigid():
    state = below_tg_limit(temperature_C=-60.0, sample_state="unpoled")
    assert np.isclose(state.x_raf, 0.4)
    assert np.isclose(state.x_maf, 0.0)
    assert np.isclose(state.roaf, 0.2)
    assert np.isclose(state.moaf, 0.0)
    assert np.isclose(state.iaf, 0.2)


def test_curve_loader_requires_traceable_figure5b_pairs(tmp_path: Path):
    path = tmp_path / "figure5b.csv"
    path.write_text(
        "figure,panel,sample_state,temperature_C,quantity,value,unit,provenance\n"
        "5,B,poled,0,x_RAF,0.16,fraction,DIGITIZED_SOURCE\n"
        "5,B,poled,0,x_MAF,0.24,fraction,DIGITIZED_SOURCE\n"
        "5,B,poled,20,x_RAF,0.12,fraction,DIGITIZED_SOURCE\n"
        "5,B,poled,20,x_MAF,0.28,fraction,DIGITIZED_SOURCE\n",
        encoding="utf-8",
    )
    states = reconstruct_curve_from_figure5b(path, sample_state="poled")
    assert len(states) == 2
    assert np.isclose(states[0].moaf, 0.04)
    assert np.isclose(states[1].moaf, 0.08)
