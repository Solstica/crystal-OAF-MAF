from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v015 import FIGURE5B_CSV, SOURCE_FIXED_FRACTIONS, V015_RULES
from pvdf_pf.calibration.source_data import load_digitized_curve_csv
from pvdf_pf.calibration.phase_fraction_state import reconstruct_curve_from_figure5b


def main() -> None:
    points = load_digitized_curve_csv(FIGURE5B_CSV)
    quantities = sorted({p.quantity for p in points})
    samples = sorted({p.sample_state.lower() for p in points})
    ready = bool(points) and {"x_RAF", "x_MAF"}.issubset(set(quantities))

    reconstructed = {}
    if ready:
        for sample in samples:
            states = reconstruct_curve_from_figure5b(
                FIGURE5B_CSV,
                sample_state=sample,
                crystal_fraction=SOURCE_FIXED_FRACTIONS["crystal"],
                iaf_fraction=SOURCE_FIXED_FRACTIONS["iaf"],
            )
            reconstructed[sample] = [s.to_dict() for s in states]

    report = {
        "model_version": "v0.1.5",
        "stage": "RAF/MAF to ROAF/MOAF phase-fraction reconstruction",
        "source_fixed_fractions": SOURCE_FIXED_FRACTIONS,
        "rules": V015_RULES,
        "figure5b_path": str(FIGURE5B_CSV.relative_to(ROOT)),
        "n_figure5b_points": len(points),
        "observed_quantities": quantities,
        "sample_states": samples,
        "ready_for_source_phase_fraction_reconstruction": ready,
        "reconstructed_states": reconstructed,
        "status": (
            "SOURCE_PHASE_FRACTIONS_READY"
            if ready
            else "WAITING_FOR_MAIN_TEXT_FIGURE_5B_DIGITIZATION"
        ),
        "guardrail": (
            "Do not infer x_RAF/x_MAF from n(T), epsilon_MOAF(T), or a guessed Tg curve. "
            "The SI explicitly states that Figure 5B values are used in Eqs. S1-S4."
        ),
    }

    outdir = ROOT / "outputs" / "v015_phase_fraction_state"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"\nsaved: {path}")


if __name__ == "__main__":
    main()
