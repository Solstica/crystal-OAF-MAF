"""v0.1.11 closed-OAF ROAF/MOAF transfer compatibility audit."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUI2022_DIR = ROOT / "data" / "literature" / "rui2022"
TEMPERATURES_C = (-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0)
SAMPLE_STATES = ("unpoled", "poled")

V0111_RULES = {
    "closed_same_state_oaf_response": "v0.1.10 epsilon_OAF,eff(T)",
    "devitrification_progress": "v0.1.6 normalized Figure-5B donor q(T); trend only",
    "mobile_oaf_local_response": "Rui-2022 Figure S2 epsilon_MOAF(T); source inverse only",
    "mixture_test": "epsilon_OAF,eff=(1-q)*epsilon_ROAF+q*epsilon_MOAF",
    "fit_new_parameters": False,
    "low_temperature_roaf_anchor_hypothesis": "epsilon_ROAF(-30C)=epsilon_OAF,eff(-30C) because q(-30C)=0 by project normalization",
    "promotion_to_tdgl": False,
}
