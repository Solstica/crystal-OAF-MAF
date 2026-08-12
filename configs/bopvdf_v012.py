"""v0.1.2 BDS / dipolar-dynamics configuration.

No numerical n(T), m_d(T), g(T), lambda(T), or mu_r(T) values are inserted here
because the accessible primary-source record currently supports the qualitative
trends and confirms that the Supporting Information contains those curves/results,
but the numerical SI dataset has not yet been transcribed into this repository.
"""

BDS_PRIMARY_SOURCE = {
    "citation": "Rui et al., Macromolecules 2022, 55, 9705-9714",
    "doi": "10.1021/acs.macromol.2c01110",
    "material": "unpoled and highly poled biaxially oriented PVDF (BOPVDF)",
    "temperature_window_C": [-30.0, 40.0],
    "supported_qualitative_constraints": [
        "active dipole concentration increases substantially from -30 to 40 C",
        "calculated dipole moment and Kirkwood-Frohlich g-factor are higher for poled BOPVDF than unpoled BOPVDF",
        "dipole-dipole interaction increases substantially with temperature",
        "rotational dipole mobility increases by more than four orders of magnitude from -30 to 40 C",
        "Supporting Information contains n(T), m_d(T), g(T), three-phase MOAF dielectric calculation, lambda(T), and mu_r(T) results",
    ],
    "phase_mapping_warning": (
        "The source interprets devitrification of the rigid amorphous fraction (RAF). "
        "RAF dynamics are not automatically identical to the oriented amorphous fraction (OAF); "
        "use these data as an amorphous/interphase dynamical constraint until an OAF-specific mapping is validated."
    ),
}

BDS_CSV_SCHEMA = {
    "required": ["frequency_Hz", "epsilon_real", "epsilon_loss"],
    "optional_grouping": ["sample_state", "temperature_C"],
    "loss_convention": "epsilon_loss is positive; internal complex convention is epsilon*=epsilon'-i*epsilon''",
}

V012_RULES = {
    "fit_model": "single_Debye_scaffold",
    "allow_direct_OAF_assignment": False,
    "allow_TDGL_time_coupling": False,
    "reason_TDGL_time_coupling_disabled": (
        "The present TDGL time step is dimensionless. A physical seconds-per-TDGL-time mapping must be calibrated before tau_s can be coupled to TDGL."
    ),
}
