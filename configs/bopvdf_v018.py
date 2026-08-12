"""v0.1.8 measured-film dielectric closure audit configuration."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUI2022_DIR = ROOT / "data" / "literature" / "rui2022"

TEMPERATURES_C = (-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0)
SAMPLE_STATES = ("unpoled", "poled")

V018_RULES = {
    "fit_new_parameters": False,
    "measured_target": "Rui 2022 main-text Figure 3A epsilon_c(T)",
    "project_prediction": "ideal parallel-capacitor closure using v0.1.6 fractions and v0.1.7 local dielectric inputs",
    "literal_si_check": "evaluate SI Eq. S4 only when eta_MOAF = x_MAF - 0.20 is non-negative",
    "on_failure": "retain residual; do not refit epsilon_MOAF or clip source fractions",
    "promotion_gate": "four-component dielectric set cannot be promoted to calibrated physical model unless film-level closure is explained",
}
