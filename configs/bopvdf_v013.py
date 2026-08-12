"""v0.1.3 source-driven OAF subpartition and BDS calibration configuration.

No plot-derived numerical data are hard-coded here.  The repository must ingest
Figure S1/S2/S3 points through the digitized CSV files before any temperature-dependent
MOAF/ROAF/IAF parameter can be promoted to a physical configuration.
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
    # Explicitly stated in Section S2 of the SI. These are source-paper model
    # assumptions, not direct local measurements.
    "epsilon_crystal": 3.0,
    "epsilon_roaf": 3.0,
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
