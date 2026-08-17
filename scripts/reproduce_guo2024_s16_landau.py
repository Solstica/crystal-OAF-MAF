"""Generate source-constrained homogeneous Landau checkpoints for Guo 2024 Fig. S16."""
from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from pvdf_pf.literature.guo2024_vector import (
    GUO2024_STRONG_ANISOTROPY_VECTOR,
    GUO2024_WEAK_ANISOTROPY_VECTOR,
    directional_anisotropy_summary,
    equilibrium_along_direction,
)


OUT = Path("outputs/literature_reproduction/guo2024_s16_landau/report.json")
T_C = 25.0
DIRECTIONS = {
    "100": np.array([1.0, 0.0, 0.0]),
    "110": np.array([1.0, 1.0, 0.0]),
    "111": np.array([1.0, 1.0, 1.0]),
}


def one_set(parameters):
    return {
        "source_parameters": asdict(parameters),
        "high_symmetry_directions": {
            name: equilibrium_along_direction(direction, T_C, parameters)
            for name, direction in DIRECTIONS.items()
        },
        "angular_anisotropy": directional_anisotropy_summary(
            T_C, parameters, n_points=4096
        ),
    }


def main() -> None:
    report = {
        "benchmark": "Guo2024 Supplementary Fig. S16 homogeneous Landau surface",
        "temperature_C": T_C,
        "source_status": {
            "coefficients": "DIRECT_TRANSCRIPTION_FROM_GUO2024_SUPPLEMENTARY_TABLES_S2_S3",
            "expanded_polynomial": "PUBLISHED_CONVENTION_CROSSCHECK_AGAINST_SU2022_EQ3",
            "figure_claim": "QUALITATIVE_ANISOTROPY_CHECKPOINT_NOT_PIXEL_EXACT_REPRODUCTION",
        },
        "strong_anisotropy": one_set(GUO2024_STRONG_ANISOTROPY_VECTOR),
        "weak_anisotropy": one_set(GUO2024_WEAK_ANISOTROPY_VECTOR),
        "interpretation_guardrail": (
            "The source provides the contracted coefficients and Fig. S16, but does not "
            "document the exact graphical radial/color normalization used for that figure. "
            "This report therefore verifies the homogeneous angular energy contrast only."
        ),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
