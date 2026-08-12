"""v0.1.5 Figure 5B source audit and ROAF/MOAF reconstruction configuration."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUI2022_DIR = ROOT / "data" / "literature" / "rui2022"
FIGURE5B_CSV = RUI2022_DIR / "figure_5b_maintext_digitized.csv"

RAW_FIGURE5B_STATE = {
    # Main-text Figure 5B: melt-recrystallized BOPVDF, xc = 0.59.
    "crystal": 0.59,
    "amorphous_total": 0.41,
}

SOURCE_FIXED_FRACTIONS = {
    # Rui et al. 2022 SI Section S2 rounded dielectric-model values.
    "crystal": 0.60,
    "iaf": 0.20,
    "total_oaf": 0.20,
}

V015_RULES = {
    "main_text_figure5b": "x_RAF + x_MAF ~= 0.41 because Figure 5B assumes xc=0.59",
    "source_equation_above_Tg": "x_RAF = eta_ROAF; x_MAF = eta_MOAF + eta_IAF",
    "source_equation_below_Tg": "x_RAF = eta_OAF + eta_IAF; x_MAF = 0",
    "si_fraction_rounding": "eta_cr=0.60, eta_IAF=0.20, eta_OAF~=0.20",
    "source_consistency_policy": (
        "audit raw Figure 5B against SI Eqs. S2-S4; do not silently force closure or "
        "clip negative eta_MOAF"
    ),
    "raf_not_oaf": True,
    "maf_not_iaf": True,
    "allow_extrapolation": False,
    "physical_morphology_activation": (
        "blocked for literal SI mapping where eta_MOAF is negative or closure fails; "
        "a separately labelled project regularization is required"
    ),
}
