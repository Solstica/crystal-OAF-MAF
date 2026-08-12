from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v015 import (
    FIGURE5B_CSV,
    RAW_FIGURE5B_STATE,
    SOURCE_FIXED_FRACTIONS,
    V015_RULES,
)
from pvdf_pf.calibration.figure5b_audit import audit_figure5b, source_anchor_checks
from pvdf_pf.calibration.source_data import load_digitized_curve_csv


def main() -> None:
    points = load_digitized_curve_csv(FIGURE5B_CSV)
    quantities = sorted({p.quantity for p in points})
    samples = sorted({p.sample_state.lower() for p in points})
    ready = bool(points) and {"x_RAF", "x_MAF"}.issubset(set(quantities))

    audit = []
    anchor_checks = {}
    if ready:
        audit = [
            p.to_dict()
            for p in audit_figure5b(
                FIGURE5B_CSV,
                sample_state="melt-recrystallized",
                raw_crystallinity=RAW_FIGURE5B_STATE["crystal"],
                si_crystal_fraction=SOURCE_FIXED_FRACTIONS["crystal"],
                si_iaf_fraction=SOURCE_FIXED_FRACTIONS["iaf"],
            )
        ]
        anchor_checks = source_anchor_checks(FIGURE5B_CSV)

    literal_valid = [p for p in audit if p["literal_si_mapping_valid"]]
    negative_moaf = [p for p in audit if not p["si_nonnegative"]]
    max_raw_closure_error = (
        max(abs(p["raw_closure_residual_vs_0p41"]) for p in audit) if audit else None
    )

    report = {
        "model_version": "v0.1.5",
        "stage": "Figure 5B source ingestion and SI consistency audit",
        "raw_figure5b_state": RAW_FIGURE5B_STATE,
        "si_fixed_fractions": SOURCE_FIXED_FRACTIONS,
        "rules": V015_RULES,
        "figure5b_path": str(FIGURE5B_CSV.relative_to(ROOT)),
        "n_figure5b_points": len(points),
        "observed_quantities": quantities,
        "sample_states": samples,
        "source_anchor_checks": anchor_checks,
        "all_source_anchors_match": bool(anchor_checks) and all(anchor_checks.values()),
        "raw_closure_max_abs_error": max_raw_closure_error,
        "n_literal_si_valid_temperatures": len(literal_valid),
        "n_negative_si_eta_moaf_temperatures": len(negative_moaf),
        "audit_points": audit,
        "ready_for_source_audit": ready,
        "ready_for_literal_si_spatial_mapping": bool(audit) and len(literal_valid) == len(audit),
        "status": (
            "SOURCE_DATA_READY_MODEL_INCONSISTENCY_FLAGGED"
            if ready and len(literal_valid) != len(audit)
            else "SOURCE_DATA_READY"
            if ready
            else "WAITING_FOR_MAIN_TEXT_FIGURE_5B_DIGITIZATION"
        ),
        "guardrail": (
            "Figure 5B is mobility-defined calorimetry with xc=0.59, while SI S2 uses rounded "
            "eta_cr=0.60 and eta_IAF=0.20. Do not clip negative eta_MOAF or force 0.40 closure "
            "without an explicitly labelled project regularization."
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
