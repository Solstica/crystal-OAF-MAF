from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pvdf_pf.literature.guo2024 import (
    GUO2024_STRONG_ANISOTROPY_AXIS,
    GUO2024_WEAK_ANISOTROPY_AXIS,
    axis_equilibrium,
    landau_axis_density,
    landau_axis_derivative,
)


def parameter_record(parameters) -> dict[str, object]:
    return {
        "label": parameters.label,
        "source_table": parameters.source_table,
        "alpha1_law": "1.412e5 * (T_C - 42) J m C^-2",
        "alpha1_prefactor_J_m_C2_K": parameters.alpha1_prefactor_J_m_C2_K,
        "curie_temperature_C": parameters.curie_temperature_C,
        "alpha11_J_m5_C4": parameters.alpha11_J_m5_C4,
        "alpha111_J_m9_C6": parameters.alpha111_J_m9_C6,
        "provenance": "DIRECT_TRANSCRIPTION_FROM_GUO2024_SUPPLEMENTARY_TABLE",
    }


def case_report(parameters, temperature_C: float) -> dict[str, object]:
    equilibrium = axis_equilibrium(temperature_C, parameters)
    p0 = equilibrium["P_abs_min_C_m2"]
    derivative_at_min = landau_axis_derivative(p0, temperature_C, parameters)
    barrier_from_min = -equilibrium["f_min_J_m3"]

    # A small table is retained so later figure reproduction can compare the
    # exact same source polynomial without fitting a new curve.
    sample_P = np.linspace(-0.10, 0.10, 81)
    sample_f = landau_axis_density(sample_P, temperature_C, parameters)
    samples = [
        {"P_C_m2": float(p), "f_J_m3": float(f)}
        for p, f in zip(sample_P, sample_f)
    ]

    return {
        "parameters": parameter_record(parameters),
        "derived_one_axis_equilibrium": {
            **equilibrium,
            "df_dP_at_positive_min": float(derivative_at_min),
            "zero_polarization_barrier_from_min_J_m3": float(barrier_from_min),
            "provenance": "PROJECT_REPRODUCTION_DERIVED_FROM_DIRECT_SOURCE_COEFFICIENTS",
        },
        "sampled_axis_curve": samples,
    }


def main() -> None:
    temperature_C = 25.0
    report = {
        "reproduction_id": "guo2024_landau_axis_stage1",
        "source": {
            "citation": "Guo et al., Nature Communications 15, 348 (2024)",
            "doi": "10.1038/s41467-023-44395-5",
            "material": "P(VDF-TrFE)",
            "source_temperature_C": temperature_C,
        },
        "scope": {
            "implemented": "one-axis Landau slice f=a1*P^2+a11*P^4+a111*P^6",
            "why_source_faithful": "all vector cross terms vanish identically on a polarization axis",
            "not_implemented": [
                "full expanded vector Landau cross-term polynomial",
                "gradient energy",
                "elastic/electrostrictive coupling",
                "electrostatics",
                "TDGL kinetics",
                "polar-spiral morphology",
            ],
            "reason_for_stopping_here": "the expanded cross-coefficient convention must be verified from the primary source/reference chain before coding it",
        },
        "strong_anisotropy": case_report(GUO2024_STRONG_ANISOTROPY_AXIS, temperature_C),
        "weak_anisotropy": case_report(GUO2024_WEAK_ANISOTROPY_AXIS, temperature_C),
    }

    outdir = ROOT / "outputs" / "literature_reproduction" / "guo2024_landau_axis"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    concise = {
        key: report[key]["derived_one_axis_equilibrium"]
        for key in ("strong_anisotropy", "weak_anisotropy")
    }
    print(json.dumps(concise, indent=2))
    print(f"saved: {path}")


if __name__ == "__main__":
    main()
