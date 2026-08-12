"""v0.1.13 same-state Figure-2 BDS spectrum calibration."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUI2022_DIR = ROOT / "data" / "literature" / "rui2022"
TEMPERATURES_C = (-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0)
SAMPLE_STATES = ("unpoled", "poled")

V0113_RULES = {
    "source_spectra": "Rui-2022 Figure 2 alpha-relaxation windows, DIGITIZED_SOURCE",
    "static_epsilon": "same-state Rui-2022 Figure 3A",
    "models_compared": ("Debye", "Cole-Cole"),
    "cole_cole_definition": "eps*=eps_inf+(eps_s-eps_inf)/(1+(i*omega*tau)^beta)",
    "alpha_CC": "1-beta",
    "loss_floor": "non-negative nuisance term only; not a material parameter",
    "crystal_fraction": 0.52,
    "crystal_epsilon": 3.0,
    "amorphous_assignment": "combined OAF+IAF only",
    "oaf_only_promotion": False,
    "right_censored_threshold": "1e7/raw_peak_frequency < 3",
    "tdgl_promotion": False,
}
