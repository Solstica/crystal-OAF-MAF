from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v0110 import RUI2022_DIR, SAMPLE_STATES, STRUCTURAL_PROXY, TEMPERATURES_C, V0110_RULES
from pvdf_pf.calibration.source_state_transfer import (
    closed_oaf_series,
    source_state_definitions,
    summarize_closed_oaf,
)


def main() -> None:
    results = {}
    for sample in SAMPLE_STATES:
        rows = closed_oaf_series(
            RUI2022_DIR,
            sample_state=sample,
            temperatures_C=TEMPERATURES_C,
        )
        results[sample] = {
            "summary": summarize_closed_oaf(rows),
            "points": [row.to_dict() for row in rows],
        }

    report = {
        "model_version": "v0.1.10",
        "stage": "source-state transfer audit and same-state closed OAF inversion",
        "rules": V0110_RULES,
        "source_states": [state.to_dict() for state in source_state_definitions()],
        "structural_proxy": STRUCTURAL_PROXY,
        "results": results,
        "status": "SAME_STATE_CLOSED_OAF_RESPONSE_READY_AS_EFFECTIVE_CONSTRAINT",
        "interpretation": [
            "The BDS target and Figure-5 TMDSC mobility donor are not the same material state: their crystallinities are about 0.52 and 0.59, respectively.",
            "The SI Figure-S2 inverse additionally rounds the donor crystallinity to 0.60 while using film epsilon(T) from the 0.52-crystalline BOPVDF target.",
            "For a same-state BOPVDF closure, Rui-2021 provides a defensible structural proxy of approximately 0.52 crystal + 0.28 OAF + 0.20 IAF.",
            "The resulting epsilon_OAF,eff(T) exactly closes the measured film response by construction and should be treated as an effective response constraint, not as an intrinsic OAF material constant.",
            "The unpoled calculation reuses the poled structural proxy only as a sensitivity comparison; the poled result is the primary source-supported use.",
        ],
    }

    outdir = ROOT / "outputs" / "v0110_source_state_transfer"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"\nsaved: {path}")


if __name__ == "__main__":
    main()
