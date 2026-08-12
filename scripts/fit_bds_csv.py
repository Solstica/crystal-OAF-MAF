from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pvdf_pf.calibration.bds import fit_debye_spectrum


REQUIRED = ("frequency_Hz", "epsilon_real", "epsilon_loss")
OPTIONAL_GROUP = ("sample_state", "temperature_C")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError("CSV has no header")
        missing = [name for name in REQUIRED if name not in reader.fieldnames]
        if missing:
            raise ValueError(f"missing required columns: {missing}")
        return list(reader)


def group_rows(rows: list[dict[str, str]]) -> dict[tuple[str, str], list[dict[str, str]]]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (row.get("sample_state", "unknown") or "unknown", row.get("temperature_C", "unknown") or "unknown")
        grouped[key].append(row)
    return dict(grouped)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fit one Debye relaxation to each BDS spectrum group. "
        "This is a spectrum-level scaffold, not a paper-specific RAF/OAF parameter extraction."
    )
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    rows = read_rows(args.csv_path)
    grouped = group_rows(rows)
    results = []
    for (sample_state, temperature_C), group in sorted(grouped.items()):
        f = np.asarray([float(row["frequency_Hz"]) for row in group])
        er = np.asarray([float(row["epsilon_real"]) for row in group])
        el = np.asarray([float(row["epsilon_loss"]) for row in group])
        fit = fit_debye_spectrum(f, er, el)
        results.append(
            {
                "sample_state": sample_state,
                "temperature_C": temperature_C,
                "fit_model": "single_Debye",
                "fit": fit.to_dict(),
                "status": "SPECTRUM_FIT_ONLY",
                "warning": (
                    "Do not identify this fitted tau directly with the Macromolecules 2022 "
                    "rotational mobility or assign the whole fitted relaxation to OAF/RAF without "
                    "the paper-specific decomposition."
                ),
            }
        )

    report = {
        "model_version": "v0.1.2",
        "input": str(args.csv_path),
        "required_columns": list(REQUIRED),
        "optional_group_columns": list(OPTIONAL_GROUP),
        "results": results,
    }
    output = args.output or args.csv_path.with_suffix(".debye_fit.json")
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"\nsaved: {output}")


if __name__ == "__main__":
    main()
