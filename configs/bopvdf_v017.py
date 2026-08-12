"""v0.1.7 four-component dielectric reconstruction configuration."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUI2022_DIR = ROOT / "data" / "literature" / "rui2022"
FIGURE5B_CSV = RUI2022_DIR / "figure_5b_maintext_digitized.csv"
FIGURE1B_CSV = RUI2022_DIR / "figure_1b_melt_eps_digitized.csv"

TEMPERATURES_C = (-30.0, 0.0, 20.0, 40.0)
SAMPLE_STATE = "poled"

SOURCE_MODEL_PERMITTIVITY = {
    "crystal": 3.0,
    "roaf": 3.0,
    "moaf": "Rui2022 Figure S2 source-model-derived epsilon_MOAF(T)",
    "iaf": (
        "Rui2021 Figure 9A/SI Sec. IV Kirkwood-Frohlich construction: infer g(T) "
        "from molten-PVDF epsilon(T), linearly extrapolate g(T), then invert for epsilon_IAF(T)"
    ),
}

V017_RULES = {
    "linear_dielectric_only": True,
    "copy_to_tdgl_eps_b": False,
    "oaf_fraction_driver": "v0.1.6 project-regularized Figure 5B devitrification progress",
    "spatial_hypotheses": ("IAF-proximal-first", "crystal-proximal-first"),
    "baseline_geometry": "source-parallel ideal laminate",
    "local_field_geometry": "45-degree periodic wavy laminate used only for morphology sensitivity",
    "iaf_source_transfer": "Kirkwood-Frohlich g(T) extrapolation; not direct epsilon(T) extrapolation",
}
