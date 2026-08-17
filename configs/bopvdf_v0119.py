"""v0.1.19 same-RMS orientation-spectrum sensitivity configuration."""
from configs.bopvdf_v0118 import (
    E_EXTERNAL_VPM,
    FROZEN_SOURCE,
    GRID,
    HUANG_SAXS_GEOMETRY_NM,
    OAF_PROFILE,
    PRONY,
    REPRESENTATIVE_STATE,
    RUI2022_DIR,
    STATIC_RELAXATION_FACTOR,
    STRUCTURE,
)

# Stay inside the v0.1.18 PROJECT_IMAGE_DERIVED_CONSTRAINT interval.
TARGET_RMS_N_ND = (0.10, 0.15, 0.20)

# Each tuple is (mode_z, relative coefficient, phase_rad).  These spatial spectra
# are project geometry hypotheses; only the target RMS normal projection inherits
# the v0.1.18 source constraint.
ORIENTATION_SPECTRA = (
    {
        "name": "single_k1",
        "modes": ((1, 1.0, 0.0),),
        "status": "PLACEHOLDER_GEOMETRY_ENSEMBLE",
    },
    {
        "name": "single_k2",
        "modes": ((2, 1.0, 0.0),),
        "status": "PLACEHOLDER_GEOMETRY_ENSEMBLE",
    },
    {
        "name": "two_mode_k1_k2",
        "modes": ((1, 1.0, 0.0), (2, 0.50, 1.0471975511965976)),
        "status": "PLACEHOLDER_GEOMETRY_ENSEMBLE",
    },
    {
        "name": "smooth_multimode_A",
        "modes": (
            (1, 1.00, 0.25),
            (2, 0.55, 1.10),
            (3, 0.30, 2.05),
            (4, 0.15, 2.70),
        ),
        "status": "PLACEHOLDER_GEOMETRY_ENSEMBLE",
    },
    {
        "name": "smooth_multimode_B",
        "modes": (
            (1, 1.00, 2.20),
            (2, 0.45, 0.40),
            (3, 0.25, 1.65),
            (4, 0.12, 2.95),
        ),
        "status": "PLACEHOLDER_GEOMETRY_ENSEMBLE",
    },
)

V0119_RULES = {
    "electrostatic_equation": "div[eps0*eps_b(r)*E + P_rel(r,t) + P_FE_frozen(r)] = 0",
    "tdgl_state": "frozen/not evolved",
    "switching_enabled": False,
    "controlled_descriptor": "sqrt(<n_ND^2>) held equal across morphology spectra",
    "controlled_descriptor_status": "0.10-0.20 lies inside v0.1.18 PROJECT_IMAGE_DERIVED_CONSTRAINT interval",
    "spatial_spectrum_status": "all harmonic mode numbers, relative amplitudes and phases are PLACEHOLDER_GEOMETRY_ENSEMBLE",
    "question": "test whether RMS normal projection alone predicts local-field RMS and upper-tail statistics",
    "dynamic_response_assignment": "Rui 2022 combined OAF+IAF generalized-Debye bank on all non-crystal cells",
    "frozen_amplitude_status": "same normalized beta/OAF sensitivity amplitudes as v0.1.16-v0.1.18",
    "oaf_profile": "uniform to isolate interface-spectrum effects",
}
