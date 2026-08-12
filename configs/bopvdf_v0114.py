"""v0.1.14 causal time-domain representation of v0.1.13 Cole-Cole fits."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUI2022_DIR = ROOT / "data" / "literature" / "rui2022"
TEMPERATURES_C = (-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0)
SAMPLE_STATES = ("unpoled", "poled")

PRONY = {
    "frequency_min_Hz": 1e-2,
    "frequency_max_Hz": 1e9,
    "n_modes": 49,
    "support_margin_decades": 1.5,
    "n_fit_frequencies": 601,
    "static_constraint_weight": 100.0,
    "max_normalized_complex_error": 5e-4,
}

V0114_RULES = {
    "input_constitutive_model": "v0.1.13 same-state combined-amorphous Cole-Cole response",
    "time_domain_representation": "positive generalized-Debye / Prony bank",
    "mode_equation": "tau_j*dP_j/dt + P_j = eps0*Delta_eps_j*E_local",
    "mode_strength_constraint": "Delta_eps_j >= 0 and sum_j Delta_eps_j = Delta_eps_amorphous",
    "mode_interpretation": "numerical quadrature modes, not separately identified molecular mechanisms",
    "frequency_validity_window_Hz": (1e-2, 1e9),
    "source_measurement_window_Hz": (1.0, 1e7),
    "response_assignment": "combined OAF+IAF amorphous response only",
    "oaf_only_promotion": False,
    "tdgl_physical_time_calibrated": False,
    "direct_tdgl_coupling": False,
}
