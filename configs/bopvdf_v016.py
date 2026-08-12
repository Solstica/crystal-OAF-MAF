"""v0.1.6 project-regularized temperature-driven OAF devitrification."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUI2022_DIR = ROOT / "data" / "literature" / "rui2022"
FIGURE5B_CSV = RUI2022_DIR / "figure_5b_maintext_digitized.csv"

TEMPERATURE_RANGE_C = (-30.0, 40.0)
STRUCTURAL_FRACTIONS = {
    "crystal": 0.60,
    "iaf": 0.20,
    "total_oaf": 0.20,
}
REFERENCE_TEMPERATURES_C = {
    "devitrification_start": -30.0,
    "devitrification_end": 40.0,
}

V016_RULES = {
    "status": "PROJECT_REGULARIZATION_HYPOTHESIS",
    "figure5b_role": "temperature-dependent mobility/devitrification progress only",
    "q_definition": "average of normalized RAF loss and normalized MAF gain from -30 to 40 C",
    "roaf_moaf_partition": "eta_ROAF=eta_OAF*(1-q); eta_MOAF=eta_OAF*q",
    "source_transfer": (
        "one Figure 5B devitrification curve is used for both unpoled and poled source states, "
        "matching the source paper's use of Figure 5B in both MOAF dielectric inversions"
    ),
    "spatial_assignment": "external mobility_score required; no physical default",
    "epsilon_moaf": "source-model-derived target; not TDGL background permittivity",
    "mu_r_to_tau": "disabled",
    "tdgl_time_coupling": "disabled",
}
