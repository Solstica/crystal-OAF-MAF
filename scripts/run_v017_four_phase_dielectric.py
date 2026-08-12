from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v017 import FIGURE1B_CSV, FIGURE5B_CSV, RUI2022_DIR, SAMPLE_STATE, TEMPERATURES_C, V017_RULES
from pvdf_pf.calibration.four_phase_dielectric import build_four_phase_permittivity_state, solve_four_phase_cell
from pvdf_pf.calibration.iaf_permittivity import fit_molten_pvdf_permittivity
from pvdf_pf.calibration.project_devitrification import build_regularized_oaf_state
from pvdf_pf.core.grid import Grid2D
from pvdf_pf.morphology.mobility_scores import crystal_proximal_score, iaf_proximal_score
from pvdf_pf.morphology.oaf_devitrification import allocate_source_oaf_subphases, source_phase_fractions
from pvdf_pf.morphology.oriented import periodic_winding_three_phase, wavy_winding_three_phase


def main() -> None:
    grid = Grid2D(nz=96, nx=96, dz=1.0, dx=1.0)
    ideal = periodic_winding_three_phase(
        grid,
        crystal_fraction=0.60,
        oaf_fraction=0.20,
        maf_fraction=0.20,
        winding_z=1,
        winding_x=0,
    )
    wavy = wavy_winding_three_phase(
        grid,
        crystal_fraction=0.60,
        oaf_fraction=0.20,
        maf_fraction=0.20,
        winding_z=1,
        winding_x=1,
        waviness_amplitude=0.06,
        waviness_mode_z=0,
        waviness_mode_x=2,
    )

    fit = fit_molten_pvdf_permittivity(FIGURE1B_CSV)
    records = []
    for T in TEMPERATURES_C:
        dev = build_regularized_oaf_state(FIGURE5B_CSV, temperature_C=T)
        frac_state = dev.as_source_fraction_state(sample_state="melt-recrystallized")
        eps_state = build_four_phase_permittivity_state(
            RUI2022_DIR,
            temperature_C=T,
            sample_state=SAMPLE_STATE,
        )

        # Source-parallel ideal laminate: the field is tangential to all layers.
        ideal_map, ideal_alloc = allocate_source_oaf_subphases(
            ideal,
            fraction_state=frac_state,
            mobility_score=iaf_proximal_score(ideal),
        )
        _, _, _, ideal_summary = solve_four_phase_cell(
            ideal_map,
            eps_state,
            grid,
            axis="x",
            E0=1.0,
        )

        sensitivity = {}
        for name, score_fn in {
            "iaf_proximal_first": iaf_proximal_score,
            "crystal_proximal_first": crystal_proximal_score,
        }.items():
            source_map, alloc = allocate_source_oaf_subphases(
                wavy,
                fraction_state=frac_state,
                mobility_score=score_fn(wavy),
            )
            _, _, _, summary = solve_four_phase_cell(
                source_map,
                eps_state,
                grid,
                axis="z",
                E0=1.0,
            )
            sensitivity[name] = {
                "allocation": alloc.__dict__,
                "realized_fractions": source_phase_fractions(source_map),
                "field_summary": summary.to_dict(),
            }

        records.append({
            "temperature_C": T,
            "devitrification": dev.to_dict(),
            "permittivity_state": eps_state.to_dict(),
            "ideal_parallel": {
                "allocation": ideal_alloc.__dict__,
                "realized_fractions": source_phase_fractions(ideal_map),
                "field_summary": ideal_summary.to_dict(),
            },
            "tilted_wavy_sensitivity": sensitivity,
        })

    report = {
        "model_version": "v0.1.7",
        "stage": "four-component source-model dielectric reconstruction",
        "rules": V017_RULES,
        "melt_figure1b_linear_fit": fit.to_dict(),
        "records": records,
        "guardrail": (
            "epsilon_crystal/ROAF/MOAF/IAF here are small-signal/source-model dielectric terms. "
            "They must not be copied directly into TDGL eps_b because orientational polarization "
            "would be double counted."
        ),
    }

    outdir = ROOT / "outputs" / "v017_four_phase_dielectric"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"\nsaved: {path}")


if __name__ == "__main__":
    main()
