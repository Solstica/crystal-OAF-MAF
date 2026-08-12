from pathlib import Path

import numpy as np

from pvdf_pf.calibration.source_data import load_digitized_curve_csv
from pvdf_pf.physics.multirelaxation import (
    RelaxationChannel,
    advance_relaxation_bank,
    total_relaxational_polarization,
)


def test_empty_digitization_template_is_valid(tmp_path: Path):
    path = tmp_path / "curve.csv"
    path.write_text(
        "figure,panel,sample_state,temperature_C,quantity,value,unit,provenance\n",
        encoding="utf-8",
    )
    assert load_digitized_curve_csv(path) == []


def test_digitized_point_requires_provenance(tmp_path: Path):
    path = tmp_path / "curve.csv"
    path.write_text(
        "figure,panel,sample_state,temperature_C,quantity,value,unit,provenance\n"
        "S2,A,poled,20,epsilon_MOAF,12.5,1,DIGITIZED_SOURCE\n",
        encoding="utf-8",
    )
    points = load_digitized_curve_csv(path)
    assert len(points) == 1
    assert points[0].quantity == "epsilon_MOAF"
    assert np.isclose(points[0].value, 12.5)


def test_relaxation_bank_keeps_channels_separate():
    shape = (3, 4)
    E = np.ones(shape) * 1e6
    states = {"moaf": np.zeros(shape), "iaf": np.zeros(shape)}
    masks = {
        "moaf": np.array([[1, 1, 0, 0], [1, 1, 0, 0], [1, 1, 0, 0]], dtype=bool),
        "iaf": np.array([[0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 1, 1]], dtype=bool),
    }
    channels = [
        RelaxationChannel("moaf", tau_s=1e-3, delta_eps=8.0),
        RelaxationChannel("iaf", tau_s=2e-3, delta_eps=5.0),
    ]
    out = advance_relaxation_bank(states, E, dt_s=1e-4, channels=channels, masks=masks)
    assert np.all(out["moaf"][~masks["moaf"]] == 0.0)
    assert np.all(out["iaf"][~masks["iaf"]] == 0.0)
    total = total_relaxational_polarization(out)
    assert np.all(total > 0.0)
