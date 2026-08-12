"""v0.1.5 source-level constants for Figure 5B versus SI rounded dielectric model."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUI2022_DIR = ROOT / "data" / "literature" / "rui2022"
FIGURE5B_REVIEWED_FILE = RUI2022_DIR / "figure_5b_maintext_digitized_READY.csv"

RAW_CALORIMETRY = {
    "crystallinity": 0.59,
    "amorphous_total": 0.41,
    "closure_tolerance": 0.012,
}

SI_DIELECTRIC_MODEL = {
    "eta_crystal": 0.60,
    "eta_iaf": 0.20,
    "eta_total_oaf": 0.20,
    "epsilon_crystal": 3.0,
    "epsilon_roaf": 3.0,
}

RULES = {
    "raw_figure5b_closure": "x_RAF + x_MAF ~= 0.41 because Figure 5B assumes x_c=0.59",
    "si_rounded_closure": "ROAF + MOAF ~= 0.20 only inside the SI rounded dielectric model",
    "raf_not_oaf": True,
    "maf_not_iaf": True,
}
