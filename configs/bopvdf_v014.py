"""v0.1.4 temperature-dependent source-state configuration."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUI2022_DIR = ROOT / "data" / "literature" / "rui2022"

TEMPERATURE_RANGE_C = (-30.0, 40.0)
SAMPLE_STATES = ("unpoled", "poled")
MOBILITY_BOUNDARIES = ("l1", "l2")

V014_RULES = {
    "ordinary_source_curve_interpolation": "piecewise linear; no extrapolation",
    "mu_r_interpolation": "piecewise linear in log10(mu_r); no extrapolation",
    "epsilon_MOAF_use": (
        "source-model-derived MOAF dielectric target; never copy directly into TDGL eps_b"
    ),
    "n_md_g_scope": "film/amorphous-model quantities; not MOAF-local",
    "lambda_mu_r_scope": "BOPVDF model quantities; not MOAF-local",
    "devitrification_index": (
        "normalized n(T) diagnostic only; not a ROAF/MOAF volume fraction"
    ),
    "roaf_moaf_fraction_status": (
        "resolved in v0.1.5 from Rui 2022 main-text Figure 5B over the source range"
    ),
    "mu_r_to_tau": "disabled until a validated constitutive relation is supplied",
    "tdgl_time_coupling": "disabled until physical seconds-per-TDGL-time is calibrated",
}
