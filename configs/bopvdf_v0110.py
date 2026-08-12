"""v0.1.10 sample-state transfer audit and closed OAF response configuration."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUI2022_DIR = ROOT / "data" / "literature" / "rui2022"
TEMPERATURES_C = (-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0)
SAMPLE_STATES = ("unpoled", "poled")

STRUCTURAL_PROXY = {
    "crystal": 0.52,
    "oaf": 0.28,
    "iaf": 0.20,
    "status": "SOURCE_ANCHORED_CLOSED_PROXY",
    "primary_scope": "poled BOPVDF",
    "warning": (
        "Rui-2021 reports crystallinity ~0.52 and OAF content ~0.28; IAF=0.20 is the closure remainder. "
        "Use as a same-state structural proxy, not as an exact phase-local volume-fraction measurement."
    ),
}

V0110_RULES = {
    "separate_bds_target_from_tmdsc_donor": True,
    "bds_target_crystallinity": 0.52,
    "tmdsc_donor_crystallinity": 0.59,
    "si_s2_rounded_crystallinity": 0.60,
    "closed_structural_proxy": STRUCTURAL_PROXY,
    "epsilon_crystal": 3.0,
    "epsilon_iaf": "Rui-2021 Kirkwood-Frohlich source-model extrapolation",
    "epsilon_oaf_effective": "inverse film-closure response; not intrinsic local OAF permittivity",
    "figure_s2_epsilon_moaf_role": "comparison/temperature-trend constraint only",
    "promotion_to_tdgl_eps_b": False,
}
