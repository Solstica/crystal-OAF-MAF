from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v013 import (
    DIGITIZED_FILES,
    MAPPING_RULES,
    PRIMARY_SOURCE,
    REQUIRED_QUANTITIES,
    SOURCE_MODEL_ASSUMPTIONS,
)
from pvdf_pf.calibration.source_data import group_by_quantity, load_digitized_curve_csv


def main() -> None:
    files = {}
    observed_quantities: set[str] = set()
    total_points = 0

    for key, path in DIGITIZED_FILES.items():
        points = load_digitized_curve_csv(path)
        grouped = group_by_quantity(points)
        observed_quantities.update(grouped)
        total_points += len(points)
        files[key] = {
            "path": str(path.relative_to(ROOT)),
            "n_points": len(points),
            "quantities": sorted(grouped),
            "required_quantities": REQUIRED_QUANTITIES[key],
        }

    required = {q for values in REQUIRED_QUANTITIES.values() for q in values}
    missing = sorted(required - observed_quantities)
    ready = len(missing) == 0

    report = {
        "model_version": "v0.1.3",
        "stage": "real-source ingestion and OAF subpartition",
        "primary_source": PRIMARY_SOURCE,
        "source_model_assumptions": SOURCE_MODEL_ASSUMPTIONS,
        "mapping_rules": MAPPING_RULES,
        "digitized_source_files": files,
        "n_digitized_points": total_points,
        "required_quantities": sorted(required),
        "missing_quantities": missing,
        "ready_for_temperature_dependent_physical_fit": ready,
        "status": "SOURCE_DATA_READY" if ready else "WAITING_FOR_DIGITIZED_SOURCE_POINTS",
        "guardrail": (
            "No numerical MOAF/ROAF/IAF dynamics will be promoted to the physical "
            "configuration until the required source curves are populated."
        ),
    }

    outdir = ROOT / "outputs" / "v013_source_ingest"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"\nsaved: {path}")


if __name__ == "__main__":
    main()
