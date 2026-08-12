from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v019 import FIGURE1B_CSV, RUI2022_DIR, SAMPLE_STATES, TEMPERATURES_C, V019_RULES
from pvdf_pf.calibration.film_closure import closure_series, summarize_closure
from pvdf_pf.calibration.iaf_permittivity import fit_molten_pvdf_permittivity


def main() -> None:
    fit = fit_molten_pvdf_permittivity(FIGURE1B_CSV)

    results = {}
    for sample in SAMPLE_STATES:
        rows = closure_series(
            RUI2022_DIR,
            sample_state=sample,
            temperatures_C=TEMPERATURES_C,
        )
        results[sample] = {
            "summary": summarize_closure(rows),
            "points": [row.to_dict() for row in rows],
        }

    report = {
        "model_version": "v0.1.9",
        "stage": "Rui-2021-to-2022 dielectric source-chain reconstruction",
        "rules": V019_RULES,
        "kirkwood_frohlich_melt_fit": fit.to_dict(),
        "results": results,
        "status": "SOURCE_FIGURE_RECONSTRUCTED_WITH_NONCLOSING_WEIGHTS",
        "interpretation": [
            "The earlier direct linear extrapolation of epsilon_IAF was incorrect; Rui 2021 extrapolates g(T) in the Kirkwood-Frohlich equation.",
            "With the corrected IAF construction, the printed Rui-2022 SI S2-S4 phase accounting still does not close consistently over the full -30..40 C window.",
            "The published Figure S2 MOAF values nearly reconstruct Figure 3A only when raw x_RAF and x_MAF are used as dielectric weights while eta_IAF=0.20 is also retained.",
            "Those source-reconstruction weights sum to about 1.21, so this numerical convention is retained only as provenance evidence and is not a physical phase-fraction field.",
            "No Figure-S2 local permittivity is promoted to TDGL eps_b or to a calibrated morphology until an admissible closed mapping is defined.",
        ],
    }

    outdir = ROOT / "outputs" / "v019_source_reconstruction"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"\nsaved: {path}")


if __name__ == "__main__":
    main()
