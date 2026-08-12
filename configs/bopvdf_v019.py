"""v0.1.9 Rui-2021/2022 source-chain reconstruction configuration."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUI2022_DIR = ROOT / "data" / "literature" / "rui2022"
FIGURE1B_CSV = RUI2022_DIR / "figure_1b_melt_eps_digitized.csv"
TEMPERATURES_C = (-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0)
SAMPLE_STATES = ("unpoled", "poled")

V019_RULES = {
    "iaf_extrapolation": (
        "Rui-2021 SI Sec. IV: infer Kirkwood-Frohlich g(T) from molten PVDF, "
        "linearly extrapolate g(T), then invert for epsilon(T)"
    ),
    "published_anchor_epsilon_am_25C": 19.5,
    "published_anchor_epsilon_am_127C": 12.1,
    "printed_si_mapping": "eta_MOAF=x_MAF-eta_IAF when non-negative",
    "source_figure_reconstruction": (
        "diagnostic only: use raw x_RAF and x_MAF as dielectric weights while retaining eta_IAF=0.20"
    ),
    "source_figure_reconstruction_is_physical_fraction_map": False,
    "promotion_to_tdgl": False,
}
