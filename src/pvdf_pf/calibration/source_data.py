"""Strict ingestion of digitized/reported literature curves.

Plot-derived points and values stated explicitly in source text remain distinguishable.
Every numeric row must retain figure/panel, sample state, unit, and provenance.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv

import numpy as np


@dataclass(frozen=True)
class DigitizedPoint:
    figure: str
    panel: str
    sample_state: str
    temperature_C: float
    quantity: str
    value: float
    unit: str
    provenance: str


REQUIRED_COLUMNS = (
    "figure",
    "panel",
    "sample_state",
    "temperature_C",
    "quantity",
    "value",
    "unit",
    "provenance",
)

ALLOWED_PROVENANCE = {
    "DIGITIZED_SOURCE",
    "DIRECT_TABULATED",
    "DIRECT_REPORTED_TEXT",
}


def load_digitized_curve_csv(path: str | Path) -> list[DigitizedPoint]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    with p.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError("digitized CSV has no header")
        missing = [name for name in REQUIRED_COLUMNS if name not in reader.fieldnames]
        if missing:
            raise ValueError(f"digitized CSV missing columns: {missing}")
        out: list[DigitizedPoint] = []
        for lineno, row in enumerate(reader, start=2):
            if all((row.get(name) or "").strip() == "" for name in REQUIRED_COLUMNS):
                continue
            if any((row.get(name) or "").strip() == "" for name in REQUIRED_COLUMNS):
                raise ValueError(f"incomplete digitized row at line {lineno}")
            try:
                temperature = float(row["temperature_C"])
                value = float(row["value"])
            except ValueError as exc:
                raise ValueError(f"non-numeric temperature/value at line {lineno}") from exc
            if not np.isfinite(temperature) or not np.isfinite(value):
                raise ValueError(f"non-finite digitized value at line {lineno}")
            provenance = row["provenance"].strip().upper()
            if provenance not in ALLOWED_PROVENANCE:
                raise ValueError(
                    f"unsupported provenance {provenance!r} at line {lineno}; "
                    f"use one of {sorted(ALLOWED_PROVENANCE)}"
                )
            out.append(
                DigitizedPoint(
                    figure=row["figure"].strip(),
                    panel=row["panel"].strip(),
                    sample_state=row["sample_state"].strip(),
                    temperature_C=temperature,
                    quantity=row["quantity"].strip(),
                    value=value,
                    unit=row["unit"].strip(),
                    provenance=provenance,
                )
            )
    return out


def group_by_quantity(points: list[DigitizedPoint]) -> dict[str, list[DigitizedPoint]]:
    grouped: dict[str, list[DigitizedPoint]] = {}
    for point in points:
        grouped.setdefault(point.quantity, []).append(point)
    for quantity in grouped:
        grouped[quantity] = sorted(
            grouped[quantity], key=lambda p: (p.sample_state, p.temperature_C)
        )
    return grouped
