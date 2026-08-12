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


def test_committed_rui2022_digitization_is_complete():
    root = Path(__file__).resolve().parents[1] / "data" / "literature" / "rui2022"
    s1 = load_digitized_curve_csv(root / "figure_s1_digitized.csv")
    s2 = load_digitized_curve_csv(root / "figure_s2_digitized.csv")
    s3 = load_digitized_curve_csv(root / "figure_s3_digitized.csv")

    assert len(s1) == 48
    assert len(s2) == 30
    assert len(s3) == 48
    assert {p.quantity for p in s1} == {"n", "m_d", "g"}
    assert {p.quantity for p in s2} == {"epsilon_MOAF"}
    assert {p.quantity for p in s3} == {"lambda", "mu_r"}
    assert all(p.provenance == "DIGITIZED_SOURCE" for p in s1 + s2 + s3)


def test_digitized_curves_reproduce_source_level_trends():
    root = Path(__file__).resolve().parents[1] / "data" / "literature" / "rui2022"
    s1 = load_digitized_curve_csv(root / "figure_s1_digitized.csv")
    s2 = load_digitized_curve_csv(root / "figure_s2_digitized.csv")
    s3 = load_digitized_curve_csv(root / "figure_s3_digitized.csv")

    for quantity in ("n", "g"):
        for temperature in range(-30, 41, 10):
            values = {
                p.sample_state: p.value
                for p in s1
                if p.quantity == quantity and p.temperature_C == temperature
            }
            assert values["poled"] > values["unpoled"]

    for temperature in range(-30, 41, 5):
        values = {
            p.sample_state: p.value
            for p in s2
            if p.quantity == "epsilon_MOAF" and p.temperature_C == temperature
        }
        assert values["poled"] > values["unpoled"]

    for sample_state in ("unpoled", "poled"):
        for panel in ("B_l1", "B_l2"):
            points = sorted(
                [
                    p
                    for p in s3
                    if p.quantity == "mu_r"
                    and p.sample_state == sample_state
                    and p.panel == panel
                ],
                key=lambda p: p.temperature_C,
            )
            assert np.log10(points[-1].value / points[0].value) > 4.0


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
