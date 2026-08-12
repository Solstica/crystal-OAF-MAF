"""Temporary promotion helper for Figure 5B source data.

This script exists only to make the source-gate transition auditable. The populated
v2 source file is copied into the canonical Figure 5B path after review; no numerical
transformation is performed here.
"""
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
src = ROOT / "data/literature/rui2022/figure_5b_maintext_digitized_v2.csv"
dst = ROOT / "data/literature/rui2022/figure_5b_maintext_digitized.csv"
if not src.exists():
    raise FileNotFoundError(src)
shutil.copyfile(src, dst)
print(f"promoted {src.relative_to(ROOT)} -> {dst.relative_to(ROOT)}")
