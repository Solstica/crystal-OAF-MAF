"""v0.1 literature-constrained BOPVDF calibration configuration.

This file does not yet define a fully physical TDGL parameter set.  It contains
macroscopic targets, phase fractions and explicit hypotheses for the first inverse
problems.  Values labelled HYPOTHESIS must not be promoted to direct measurements.
"""

from pvdf_pf.calibration.targets import (
    BOPVDF_2021_TARGETS,
    MACROMOLECULES_2022_ASSUMPTIONS,
    YANG_2015_TARGETS,
)

GRID = {"nz": 200, "nx": 32, "dz": 1.0, "dx": 1.0}

PHASE_FRACTIONS = {
    "crystal": BOPVDF_2021_TARGETS["crystal_fraction"]["value"],
    "oaf": BOPVDF_2021_TARGETS["oaf_fraction"]["value"],
    "maf": BOPVDF_2021_TARGETS["maf_fraction"]["value"],
}

MORPHOLOGY = {
    "period": 100,
    "offset": 0,
}

DIELECTRIC_CALIBRATION = {
    "target_eps_eff": BOPVDF_2021_TARGETS["epsilon_r_low_field"]["value"],
    # The crystal value below is a SOURCE_MODEL_ASSUMPTION from the 2022
    # Macromolecules supporting information, not a direct local measurement.
    "crystal_eps_r": MACROMOLECULES_2022_ASSUMPTIONS["crystal_epsilon_r"]["value"],
    # Yang 2015 reports 21-22 for the amorphous phase.  v0.1 deliberately uses
    # its midpoint only as a sensitivity hypothesis for the MAF; OAF/MAF were
    # not independently measured in that extraction.
    "maf_eps_hypothesis": YANG_2015_TARGETS["amorphous_phase_epsilon_r"]["value"],
    "maf_eps_scan": [8.0, 12.0, 16.0, 21.5, 26.0],
    "oaf_fit_bracket": [1.01, 500.0],
    # The phase map varies along z.  x therefore corresponds to a field parallel
    # to idealized lamellae; z corresponds to a field normal to them.
    "axes": ["x", "z"],
}

POLARIZATION_CALIBRATION = {
    "film_ps_Cpm2": BOPVDF_2021_TARGETS["spontaneous_polarization_Cpm2"]["value"],
    "beta_crystal_ps_upper_bound_Cpm2": BOPVDF_2021_TARGETS["beta_crystal_Ps_upper_bound_Cpm2"]["value"],
    "maf_ps_first_bound_Cpm2": 0.0,
}
