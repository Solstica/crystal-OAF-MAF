"""v0.1.3 source-driven OAF subpartition and BDS calibration configuration.

Plot-derived numerical data are ingested only through the digitized CSV files.
The values in SOURCE_MODEL_ASSUMPTIONS are the assumptions/partition used by
Rui et al. 2022 SI; they are not universal local material constants and they do
not overwrite the earlier v0.1 project baseline fractions.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUI2022_DIR = ROOT / "data" / "literature" / "rui2022"

PRIMARY_SOURCE = {
    "citation": "Rui et al., Macromolecules 2022, 55, 9705-9714",
    "doi": "10.1021/acs.macromol.2c01110",
    "supporting_information": "ma2c01110_si_001.pdf",
    "temperature_window_C": [-30.0, 40.0],
}

SOURCE_MODEL_ASSUMPTIONS = {
    # Explicitly stated in Sections S2-S3 of the SI. These are source-paper
    # model assumptions/partitions, not direct phase-resolved measurements.
    "epsilon_crystal": 3.0,
    "epsilon_roaf": 3.0,
    "eta_crystal": 0.6,
    "eta_iaf_constant": 0.2,
    "eta_oaf_approx": 0.2,
    "eta_oaf_relation": "eta_ROAF(T) + eta_MOAF(T)",
    "below_Tg": {
        "x_RAF": "eta_OAF(T) + eta_IAF(T)",
        "x_MAF": 0.0,
    },
    "above_Tg": {
        "x_RAF": "eta_ROAF(T)",
        "x_MAF": "eta_MOAF(T) + eta_IAF(T)",
    },
}

DIGITIZED_FILES = {
    "figure_s1": RUI2022_DIR / "figure_s1_digitized.csv",
    "figure_s2": RUI2022_DIR / "figure_s2_digitized.csv",
    "figure_s3": RUI2022_DIR / "figure_s3_digitized.csv",
}

REQUIRED_QUANTITIES = {
    "figure_s1": ["n", "m_d", "g"],
    "figure_s2": ["epsilon_MOAF"],
    "figure_s3": ["lambda", "mu_r"],
}

MAPPING_RULES = {
    "project_total_OAF": "ROAF + MOAF",
    "project_MAF_to_source_IAF": "provisional mapping only; validate before physical use",
    "allow_interpolation_of_source_points": False,
    "allow_unlabeled_plot_values": False,
    "allow_direct_TDGL_time_coupling": False,
}
