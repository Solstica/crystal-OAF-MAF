"""v0.1.16 frozen-ferroelectric-source orientation audit configuration.

The Huang 2021 SAXS layer thicknesses are used only to set geometric volume
fractions.  The WAXD crystallinity 0.52 is a weight fraction and is retained as a
source datum rather than substituted as voxel occupancy.
"""
from configs.bopvdf_v0115 import PRONY, REPRESENTATIVE_STATE, RUI2022_DIR

HUANG_SAXS_GEOMETRY_NM = {
    "lamellar_period": 11.8,
    "beta_crystal": 5.78,
    "oaf_total": 3.02,
    "oaf_each_side": 1.51,
    "iaf": 3.00,
}

L = HUANG_SAXS_GEOMETRY_NM["lamellar_period"]
STRUCTURE = {
    "crystal_fraction_volume": HUANG_SAXS_GEOMETRY_NM["beta_crystal"] / L,
    "oaf_fraction_volume": HUANG_SAXS_GEOMETRY_NM["oaf_total"] / L,
    "iaf_fraction_volume": HUANG_SAXS_GEOMETRY_NM["iaf"] / L,
    "source_beta_weight_fraction_waxd": 0.52,
    "source_min_oaf_weight_fraction_model": 0.25,
    "epsilon_crystal_background": 3.0,
}

# Square normalized cell: the electrostatic audit depends on layer orientation and
# volume fractions; no gradient-energy length scale is used in this version.
GRID = {
    "nz": 48,
    "nx": 48,
    "dz": 1.0,
    "dx": 1.0,
}

# Deliberately normalized source magnitudes.  Absolute beta/OAF remanent
# polarization is underidentified by the available bulk Pr0 data, so these values
# are not presented as measurements.  The linear electrostatic response can be
# rescaled once a phase allocation is calibrated.
FROZEN_SOURCE = {
    "p_beta_C_m2": 0.010,
    "p_oaf_mean_C_m2": 0.006,
    "oaf_decay_length_fraction_of_half_oaf": 0.35,
}

ORIENTATIONS = (
    {"name": "huang_aligned_tangent", "winding_z": 0, "winding_x": 1},
    {"name": "tilted_45deg_control", "winding_z": 1, "winding_x": 1},
    {"name": "normal_to_layers_control", "winding_z": 1, "winding_x": 0},
)

OAF_PROFILES = ("uniform", "interface_decay")
STATIC_RELAXATION_FACTOR = 100.0
E_EXTERNAL_VPM = 0.0

V0116_RULES = {
    "electrostatic_equation": "div[eps0*eps_b(r)*E + P_rel(r,t) + P_FE_frozen(r)] = 0",
    "tdgl_state": "frozen/not evolved",
    "switching_enabled": False,
    "huang_geometry_transfer": "SAXS layer thickness ratios set voxel volume fractions",
    "weight_volume_distinction": "WAXD beta fraction 0.52 is retained as a weight-fraction datum and is not used as voxel occupancy",
    "polarization_orientation": "scalar P_FE,z follows the film-normal poling direction",
    "huang_aligned_interface_relation": "SAXS lamellar normal is along MD, so film-normal P is tangent to ideal flat interfaces",
    "oaf_profile_status": "uniform and interface-decay profiles are sensitivity hypotheses with identical OAF mean polarization",
    "frozen_amplitude_status": "normalized sensitivity values; not a phase-resolved remanent-polarization measurement",
    "dynamic_response_assignment": "Rui 2022 combined OAF+IAF generalized-Debye bank on all non-crystal cells",
    "cross_source_transfer": "Huang 2021 geometry is combined with the already calibrated Rui 2022 BDS constitutive bank only for an orientation/electrostatic sensitivity audit",
}
