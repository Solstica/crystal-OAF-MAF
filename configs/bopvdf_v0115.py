"""v0.1.15 self-consistent local-field coupling configuration."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUI2022_DIR = ROOT / "data" / "literature" / "rui2022"

REPRESENTATIVE_STATE = {
    "sample_state": "poled",
    "temperature_C": 0.0,
}

STRUCTURE = {
    "crystal_fraction": 0.52,
    "oaf_fraction": 0.28,
    "maf_fraction": 0.20,
    "epsilon_crystal_background": 3.0,
}

GRID = {
    "nz": 50,
    "nx": 50,
    "dz": 1.0,
    "dx": 1.0,
}

E_EXTERNAL_VPM = 1.0e6
STATIC_RELAXATION_FACTOR = 100.0
DYNAMIC_REFERENCE_FREQUENCIES_HZ = (1.0, 1.0e3, 1.0e6)

PRONY = {
    "frequency_min_Hz": 1e-2,
    "frequency_max_Hz": 1e9,
    "n_modes": 49,
    "support_margin_decades": 1.5,
    "n_fit_frequencies": 601,
    "static_constraint_weight": 100.0,
}

V0115_RULES = {
    "electrostatic_equation": "div[eps0*eps_b(r)*E(r,t)+P_rel(r,t)] = 0",
    "applied_field_treatment": "periodic corrector includes dielectric heterogeneity under uniform macroscopic Ez",
    "dynamic_response_assignment": "same generalized-Debye bank on all non-crystal cells (combined OAF+IAF response)",
    "crystal_background_epsilon": 3.0,
    "amorphous_background_epsilon": "v0.1.13 fitted epsilon_infinity_amorphous(T)",
    "time_step_coupling": "algebraically self-consistent zero-order-hold elimination; one linear Poisson solve per step",
    "tdgl_state": "frozen/not coupled",
    "oaf_iaf_separated_constitutive_laws": False,
    "physical_tdgl_time_calibrated": False,
}
