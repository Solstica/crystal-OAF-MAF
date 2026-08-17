"""v0.1.17 lamellar-waviness sensitivity configuration.

The Huang 2021 SAXS orientation is retained on average: the ideal lamellar normal
is along MD while the frozen polarization is along the film normal.  Sinusoidal
waviness locally tilts the beta/OAF/IAF interfaces and therefore gives their
normal a film-normal component without rotating the whole laminate.

No experimental azimuthal-width or interface-waviness distribution is available
in the current source set, so every amplitude below is a sensitivity parameter.
"""
from configs.bopvdf_v0116 import (
    E_EXTERNAL_VPM,
    FROZEN_SOURCE,
    GRID,
    HUANG_SAXS_GEOMETRY_NM,
    PRONY,
    REPRESENTATIVE_STATE,
    RUI2022_DIR,
    STATIC_RELAXATION_FACTOR,
    STRUCTURE,
)

# Waviness amplitude is expressed as a fraction of one lamellar period in the
# periodic morphology coordinate u = x/Lx + A sin(2*pi*z/Lz).
WAVINESS_AMPLITUDES = (0.0, 0.01, 0.02, 0.05, 0.10, 0.15)
WAVINESS_MODE_Z = 1
WAVINESS_MODE_X = 0

# The first audit keeps the phase-resolved frozen source uniform inside each OAF
# layer.  This isolates geometric tilt/waviness from the separate OAF-profile
# hypothesis tested in v0.1.16.
OAF_PROFILE = "uniform"

V0117_RULES = {
    "electrostatic_equation": "div[eps0*eps_b(r)*E + P_rel(r,t) + P_FE_frozen(r)] = 0",
    "tdgl_state": "frozen/not evolved",
    "switching_enabled": False,
    "mean_orientation": "Huang-aligned: lamellar normal along MD, polarization along film normal",
    "waviness_geometry": "u = x/Lx + A*sin(2*pi*z/Lz); A is a fraction of one lamellar period",
    "waviness_status": "uncalibrated sensitivity parameter; current Huang/Rui sources do not report a numerical azimuthal FWHM or interface-waviness distribution",
    "physical_test": "local interface tilt gives n_z != 0 and therefore sigma_b = Delta(P_FE)*n_z at beta/OAF/IAF boundaries",
    "source_amplitude_status": "same normalized beta/OAF sensitivity amplitudes as v0.1.16",
    "dynamic_response_assignment": "Rui 2022 combined OAF+IAF generalized-Debye bank on all non-crystal cells",
    "oaf_profile": "uniform only in v0.1.17 to isolate morphology waviness",
}
