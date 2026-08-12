"""v0.1.12 same-state BOPVDF BDS calibration rules."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUI2022_DIR = ROOT / "data" / "literature" / "rui2022"
TEMPERATURES_C = (-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0)
SAMPLE_STATES = ("unpoled", "poled")

ETA_CRYSTAL = 0.52
EPSILON_CRYSTAL = 3.0

REPORTED_ARRHENIUS_EA_KJ_MOL = {
    "unpoled": {"mean": 57.8, "uncertainty": 1.5},
    "poled": {"mean": 68.3, "uncertainty": 3.6},
}

V0112_RULES = {
    "film_static_permittivity": "Rui-2022 Figure 3A digitized epsilon_c(T)",
    "alpha_peak_frequency": "Rui-2022 Figure 3B digitized ln(f_d)",
    "crystal_fraction": "eta_cr=0.52; direct article text for unpoled and poled BOPVDF",
    "crystal_permittivity": "epsilon_cr=3; direct article Eq. 20 treatment",
    "amorphous_static_response": "epsilon_am=(epsilon_film-eta_cr*epsilon_cr)/(1-eta_cr)",
    "dynamic_assignment": "combined OAF+IAF amorphous alpha relaxation only",
    "oaf_only_tau": False,
    "oaf_only_complex_permittivity": False,
    "fast_limit_amorphous_permittivity": "UNRESOLVED_IN_V0.1.12",
    "debye_bridge": "available only when an independent epsilon_amorphous_fast is supplied",
    "promotion_to_tdgl": False,
}
