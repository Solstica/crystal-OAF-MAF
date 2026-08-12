from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v01 import (
    DIELECTRIC_CALIBRATION,
    GRID,
    MORPHOLOGY,
    PHASE_FRACTIONS,
    POLARIZATION_CALIBRATION,
)
from pvdf_pf.calibration.dielectric import infer_unknown_phase_permittivity
from pvdf_pf.calibration.polarization import oaf_polarization_lower_bound
from pvdf_pf.calibration.targets import validate_fraction_closure
from pvdf_pf.core.grid import Grid2D
from pvdf_pf.morphology.three_phase import lamellar_three_phase_from_fractions, phase_fractions


def main() -> None:
    validate_fraction_closure()
    grid = Grid2D(**GRID)
    phase_map = lamellar_three_phase_from_fractions(
        grid,
        crystal_fraction=PHASE_FRACTIONS["crystal"],
        oaf_fraction=PHASE_FRACTIONS["oaf"],
        maf_fraction=PHASE_FRACTIONS["maf"],
        **MORPHOLOGY,
    )

    realized = phase_fractions(phase_map)
    p_oaf_min = oaf_polarization_lower_bound(
        film_ps_Cpm2=POLARIZATION_CALIBRATION["film_ps_Cpm2"],
        crystal_fraction=realized["crystal"],
        oaf_fraction=realized["oaf"],
        crystal_ps_upper_bound_Cpm2=POLARIZATION_CALIBRATION["beta_crystal_ps_upper_bound_Cpm2"],
        maf_ps_Cpm2=POLARIZATION_CALIBRATION["maf_ps_first_bound_Cpm2"],
        maf_fraction=realized["maf"],
    )

    dielectric_results = []
    for maf_eps in DIELECTRIC_CALIBRATION["maf_eps_scan"]:
        for axis in DIELECTRIC_CALIBRATION["axes"]:
            row = {
                "axis": axis,
                "crystal_eps_r": DIELECTRIC_CALIBRATION["crystal_eps_r"],
                "maf_eps_r_hypothesis": maf_eps,
                "target_eps_eff": DIELECTRIC_CALIBRATION["target_eps_eff"],
            }
            try:
                row["oaf_eps_r_inferred"] = infer_unknown_phase_permittivity(
                    phase_map,
                    grid,
                    target_eps_eff=DIELECTRIC_CALIBRATION["target_eps_eff"],
                    known_phase_eps={
                        "crystal": DIELECTRIC_CALIBRATION["crystal_eps_r"],
                        "maf": maf_eps,
                    },
                    unknown_phase="oaf",
                    axis=axis,
                    bracket=tuple(DIELECTRIC_CALIBRATION["oaf_fit_bracket"]),
                )
                row["status"] = "FIT"
            except ValueError as exc:
                row["oaf_eps_r_inferred"] = None
                row["status"] = "TARGET_UNREACHABLE"
                row["reason"] = str(exc)
            dielectric_results.append(row)

    report = {
        "model_version": "v0.1",
        "purpose": "literature-constrained parameter bounds; not a fully calibrated TDGL model",
        "phase_fraction_target": PHASE_FRACTIONS,
        "phase_fraction_realized": realized,
        "polarization_constraint": {
            "film_Ps_Cpm2": POLARIZATION_CALIBRATION["film_ps_Cpm2"],
            "beta_crystal_Ps_upper_bound_Cpm2": POLARIZATION_CALIBRATION["beta_crystal_ps_upper_bound_Cpm2"],
            "MAF_Ps_assumed_Cpm2": POLARIZATION_CALIBRATION["maf_ps_first_bound_Cpm2"],
            "OAF_local_Ps_lower_bound_Cpm2": p_oaf_min,
            "interpretation": "lower bound obtained by assigning the crystal its literature upper bound; lower crystal Ps requires higher OAF Ps",
        },
        "dielectric_inverse_problem": dielectric_results,
        "warnings": [
            "The 21-22 amorphous-phase permittivity reported by Yang 2015 is not treated as a direct MAF measurement.",
            "The crystal eps_r=3 value is a source-paper model assumption from Macromolecules 2022 SI, not a direct local measurement.",
            "Inferred OAF permittivity is morphology- and orientation-dependent and is not a measured material constant.",
            "Measured small-signal phase permittivity must not be copied into TDGL eps_b without separating explicit dipolar polarization, or dipolar response will be double counted.",
        ],
    }

    outdir = ROOT / "outputs" / "v01_calibration"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "calibration_report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(json.dumps(report, indent=2))
    print(f"\nsaved: {path}")


if __name__ == "__main__":
    main()
