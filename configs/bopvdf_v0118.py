"""v0.1.18 Huang-SAXS orientation-constrained geometry scan."""
from configs.bopvdf_v0117 import (
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
    WAVINESS_MODE_X,
    WAVINESS_MODE_Z,
)

# Control variable for v0.1.18.  These values describe the RMS film-normal
# projection of the local lamellar/interface normal, not a directly reported
# Huang et al. observable.
SOURCE_CONSTRAINED_RMS_NZ = (0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20)
SENSITIVITY_BOUND_RMS_NZ = (0.25, 0.30)
RMS_NZ_SCAN = (0.0,) + SOURCE_CONSTRAINED_RMS_NZ + SENSITIVITY_BOUND_RMS_NZ

ORIENTATION_PROVENANCE = {
    "source": "Huang 2021 Supplementary Fig. 13a rasterized 2D edge-on SAXS pattern",
    "directly_reported_azimuthal_profile": False,
    "directly_reported_azimuthal_fwhm": False,
    "image_estimated_fwhm_deg": (9.0, 20.0),
    "image_estimated_gaussian_sigma_deg": (3.8, 8.5),
    "source_constrained_rms_nz_interval": (0.05, 0.20),
    "source_constrained_status": "PROJECT_IMAGE_DERIVED_CONSTRAINT",
    "upper_sensitivity_status": "PLACEHOLDER_SENSITIVITY_BOUND",
    "replacement_rule": "replace with raw/tabulated I(chi), Herman factor, numerical azimuthal FWHM, or resolved real-space normal distribution when available",
}

V0118_RULES = {
    "electrostatic_equation": "div[eps0*eps_b(r)*E + P_rel(r,t) + P_FE_frozen(r)] = 0",
    "tdgl_state": "frozen/not evolved",
    "switching_enabled": False,
    "geometry_control_variable": "sqrt(<n_ND^2>) = RMS film-normal projection of local interface normal",
    "amplitude_role": "internal sinusoidal amplitude solved numerically from target RMS n_ND; not the reported scan variable",
    "source_constraint_status": "0.05-0.20 is image-derived from rasterized Huang Supplementary Fig. 13a and is not a direct reported measurement",
    "sensitivity_bound_status": "0.25 and 0.30 are explicit upper sensitivity bounds",
    "frozen_amplitude_status": "same normalized beta/OAF sensitivity amplitudes as v0.1.16-v0.1.17",
    "dynamic_response_assignment": "Rui 2022 combined OAF+IAF generalized-Debye bank on all non-crystal cells",
    "oaf_profile": "uniform to keep morphology orientation spread isolated",
}
