from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v0111 import RUI2022_DIR, SAMPLE_STATES, TEMPERATURES_C, V0111_RULES
from pvdf_pf.calibration.oaf_compatibility import compatibility_series, summarize_compatibility


def main() -> None:
    results = {}
    for sample in SAMPLE_STATES:
        rows = compatibility_series(
            RUI2022_DIR,
            sample_state=sample,
            temperatures_C=TEMPERATURES_C,
        )
        results[sample] = {
            "summary": summarize_compatibility(rows),
            "points": [row.to_dict() for row in rows],
        }

    report = {
        "model_version": "v0.1.11",
        "stage": "closed-OAF ROAF/MOAF transfer compatibility audit",
        "rules": V0111_RULES,
        "results": results,
        "status": "TRANSFER_COMPATIBILITY_AUDITED_NO_FORCED_FIT",
        "interpretation": [
            "The Figure-5B normalized devitrification progress and Figure-S2 epsilon_MOAF are transferred from a source construction that is not the same structural state as the BOPVDF film target.",
            "At q<1, the audit reports the epsilon_ROAF(T) required to satisfy the closed v0.1.10 OAF response exactly; strong temperature drift means a constant rigid-OAF response is not supported by the combined transferred quantities.",
            "At q=1, epsilon_ROAF drops out of the mixture equation. Any remaining epsilon_OAF,eff-epsilon_MOAF gap is therefore an irreducible endpoint incompatibility of the transferred q and Figure-S2 epsilon_MOAF under the closed same-state proxy.",
            "A constant-low-T ROAF anchor is evaluated only as a diagnostic. Its missing response is not introduced as a fitted material term.",
            "The next admissible calibration must either determine BOPVDF-specific OAF mobility fractions/relaxation response from same-state BDS data or introduce an independently justified collective/interfacial response; the present audit does neither automatically.",
        ],
    }

    outdir = ROOT / "outputs" / "v0111_oaf_compatibility"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"\nsaved: {path}")


if __name__ == "__main__":
    main()
