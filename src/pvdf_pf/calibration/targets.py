"""Literature targets used by v0.1 calibration.

The entries intentionally distinguish measured quantities, model-derived phase
fractions, and assumptions introduced by the source papers.  No local OAF or MAF
permittivity is labelled as directly measured unless the source actually reports it.
"""

BOPVDF_2021_TARGETS = {
    "material": "highly poled biaxially oriented PVDF (BOPVDF), pure beta crystals",
    "temperature_K": 298.15,
    "frequency_Hz": 10.0,
    "field_for_loop_Vpm": 300e6,
    "crystal_fraction": {
        "value": 0.52,
        "status": "DIRECT",
        "method": "WAXD crystallinity",
        "source": "Huang et al., Nature Communications 12, 675 (2021), doi:10.1038/s41467-020-20662-7",
    },
    "oaf_fraction": {
        "value": 0.28,
        "status": "INFERRED_FROM_HYSTERESIS",
        "method": "three-phase analysis of poled BOPVDF hysteresis",
        "source": "Rui et al., J. Mater. Chem. C 9, 894 (2021), doi:10.1039/D0TC04632A",
    },
    "maf_fraction": {
        "value": 0.20,
        "status": "INFERRED_BY_CLOSURE",
        "method": "1 - f_crystal - f_OAF",
        "source": "derived from the two entries above",
    },
    "epsilon_r_low_field": {
        "value": 19.5,
        "status": "DIRECT",
        "method": "BDS, real relative permittivity at 25 C and 10 Hz",
        "source": "Huang et al., Nature Communications 12, 675 (2021)",
    },
    "epsilon_r_fresh_reference": {
        "value": 11.5,
        "status": "DIRECT",
        "method": "BDS, fresh BOPVDF at 25 C and 10 Hz",
        "source": "Huang et al., Nature Communications 12, 675 (2021)",
    },
    "dynamic_permittivity_high_field": {
        "value": 22.9,
        "status": "DIRECT",
        "method": "slope dD/d(eps0 E) of poled BOPVDF at high field",
        "source": "Huang et al., Nature Communications 12, 675 (2021)",
    },
    "spontaneous_polarization_Cpm2": {
        "value": 0.140,
        "status": "DIRECT",
        "method": "D-E loop at 300 MV/m, 10 Hz",
        "source": "Huang et al., Nature Communications 12, 675 (2021); Rui et al., J. Mater. Chem. C 9, 894 (2021)",
    },
    "coercive_field_Vpm": {
        "value": 88e6,
        "status": "DIRECT",
        "method": "D-E loop",
        "source": "Huang et al., Nature Communications 12, 675 (2021)",
    },
    "beta_crystal_Ps_upper_bound_Cpm2": {
        "value": 0.188,
        "status": "ATOMISTIC_REFERENCE",
        "method": "DFT theoretical limit cited and used by Huang et al.",
        "source": "Huang et al., Nature Communications 12, 675 (2021)",
    },
}


YANG_2015_TARGETS = {
    "material": "commercial BOPVDF",
    "temperature_K": 298.15,
    "amorphous_phase_epsilon_r": {
        "value": 21.5,
        "range": [21.0, 22.0],
        "status": "DIRECT_REPORTED_RANGE",
        "method": "extracted amorphous-phase dielectric constant",
        "source": "Yang et al., ACS Appl. Mater. Interfaces 7, 19894-19905 (2015), doi:10.1021/acsami.5b02944",
        "note": "Do not silently equate this value with a pure MAF local permittivity; the 2015 analysis motivated a three-phase lamellar crystal/oriented interphase/amorphous model.",
    },
}


MACROMOLECULES_2022_ASSUMPTIONS = {
    "material": "unpoled and highly poled BOPVDF",
    "crystal_epsilon_r": {
        "value": 3.0,
        "status": "SOURCE_MODEL_ASSUMPTION",
        "source": "Rui & Allahyarov et al., Macromolecules 55, 9705-9714 (2022), doi:10.1021/acs.macromol.2c01110, Supporting Information S2",
    },
    "rigid_oaf_epsilon_r": {
        "value": 3.0,
        "status": "SOURCE_MODEL_ASSUMPTION",
        "source": "same Supporting Information S2; epsilon_ROAF assumed equal to epsilon_crystal",
    },
    "iaf_epsilon_r": {
        "value": None,
        "status": "SOURCE_MODEL_DERIVED",
        "source": "same Supporting Information S2; obtained by extrapolation from melt permittivity",
    },
    "mobile_oaf_epsilon_r": {
        "value": None,
        "status": "SOURCE_MODEL_DERIVED",
        "source": "same Supporting Information S2; calculated with a temperature-dependent three-phase model",
    },
    "dynamic_variables": {
        "names": ["active_dipole_concentration", "dipole_moment", "Kirkwood_Frohlich_g", "rotational_dipole_mobility"],
        "status": "DIRECT_MODEL_OUTPUT_FROM_BDS",
        "source": "Rui & Allahyarov et al., Macromolecules 55, 9705-9714 (2022)",
    },
}


def validate_fraction_closure(targets: dict = BOPVDF_2021_TARGETS, atol: float = 1e-12) -> float:
    total = sum(targets[key]["value"] for key in ("crystal_fraction", "oaf_fraction", "maf_fraction"))
    if abs(total - 1.0) > atol:
        raise ValueError(f"phase fractions do not close: sum={total}")
    return total
