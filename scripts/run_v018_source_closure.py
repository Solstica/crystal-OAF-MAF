from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v018 import RUI2022_DIR, SAMPLE_STATES, TEMPERATURES_C, V018_RULES
from pvdf_pf.calibration.film_closure import closure_series, summarize_closure


def main() -> None:
    results = {}
    for sample in SAMPLE_STATES:
        rows = closure_series(
            RUI2022_DIR,
            sample_state=sample,
            temperatures_C=TEMPERATURES_C,
        )
        results[sample] = {
            "summary": summarize_closure(rows),
            "points": [r.to_dict() for r in rows],
        }

    max_rel = max(
        results[s]["summary"]["max_abs_relative_error_project"] for s in SAMPLE_STATES
    )
    status = (
        "SOURCE_LOCAL_DIELECTRIC_SET_CLOSES_TO_FILM_DATA"
        if max_rel <= 0.05
        else "SOURCE_LOCAL_DIELECTRIC_SET_NOT_CLOSED_TO_FILM_DATA"
    )

    report = {
        "model_version": "v0.1.8",
        "stage": "measured-film dielectric closure audit",
        "rules": V018_RULES,
        "status": status,
        "results": results,
        "interpretation_guardrail": (
            "A closure failure is not repaired by refitting source epsilon_MOAF or by silently "
            "renormalizing Figure 5B. It means the present combination of source assumptions, "
            "source-derived local permittivities, project regularized OAF fractions, and ideal "
            "parallel accounting is not yet a quantitatively closed film model."
        ),
    }

    outdir = ROOT / "outputs" / "v018_source_closure"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"\nsaved: {path}")


if __name__ == "__main__":
    main()
