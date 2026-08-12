"""v0.1.5 ROAF/MOAF phase-fraction reconstruction configuration."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUI2022_DIR = ROOT / "data" / "literature" / "rui2022"
FIGURE5B_CSV = RUI2022_DIR / "figure_5b_maintext_digitized.csv"

SOURCE_FIXED_FRACTIONS = {
    # Rui et al. 2022 SI Section S2 source-model values.
    "crystal": 0.6,
    "iaf": 0.2,
    "total_oaf": 0.2,
}

V015_RULES = {
    "source_equation_above_Tg": "x_RAF = eta_ROAF; x_MAF = eta_MOAF + eta_IAF",
    "source_equation_below_Tg": "x_RAF = eta_OAF + eta_IAF; x_MAF = 0",
    "figure5b_required": True,
    "allow_figure5b_inference_from_nT": False,
    "allow_figure5b_inference_from_epsilon_MOAF": False,
    "allow_extrapolation": False,
    "source_fraction_set_scope": (
        "Rui 2022 SI source model only; do not merge numerically with the separate "
        "0.52/0.28/0.20 BOPVDF parameter set"
    ),
    "physical_morphology_activation": (
        "disabled until Figure 5B x_RAF(T)/x_MAF(T) points are traceably digitized"
    ),
}
