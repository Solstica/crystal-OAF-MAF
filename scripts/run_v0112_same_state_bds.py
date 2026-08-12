from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v0112 import (
    EPSILON_CRYSTAL,
    ETA_CRYSTAL,
    REPORTED_ARRHENIUS_EA_KJ_MOL,
    RUI2022_DIR,
    SAMPLE_STATES,
    TEMPERATURES_C,
    V0112_RULES,
)
from pvdf_pf.calibration.same_state_bds import (
    fit_high_temperature_arrhenius,
    same_state_bds_series,
)


def main() -> None:
    results = {}
    for sample in SAMPLE_STATES:
        rows = same_state_bds_series(
            RUI2022_DIR,
            sample_state=sample,
            temperatures_C=TEMPERATURES_C,
            eta_crystal=ETA_CRYSTAL,
            epsilon_crystal=EPSILON_CRYSTAL,
        )
        arr = fit_high_temperature_arrhenius(rows)
        reported = REPORTED_ARRHENIUS_EA_KJ_MOL[sample]
        rel_error = abs(arr.activation_energy_kj_mol - reported["mean"]) / reported["mean"]
        results[sample] = {
            "arrhenius_digitization_check": {
                **arr.to_dict(),
                "reported_activation_energy_kj_mol": reported["mean"],
                "reported_uncertainty_kj_mol": reported["uncertainty"],
                "relative_difference_from_reported_mean": rel_error,
            },
            "points": [row.to_dict() for row in rows],
        }

    report = {
        "model_version": "v0.1.12",
        "stage": "same-state BOPVDF broadband-dielectric relaxation bridge",
        "rules": V0112_RULES,
        "results": results,
        "status": "SAME_STATE_COMBINED_AMORPHOUS_BDS_ESTABLISHED",
        "interpretation": [
            "Figure 3A static permittivity and Figure 3B alpha-relaxation peak frequency come from the same unpoled/poled BOPVDF film states.",
            "Equation 20 with eta_cr=0.52 and epsilon_cr=3 yields the combined amorphous static permittivity epsilon_am(T).",
            "The article explicitly treats OAF and IAF together for this BDS inversion, so v0.1.12 does not assign the measured alpha-relaxation time to OAF alone.",
            "Digitized positive-temperature Figure 3B points are checked against the article-reported Arrhenius activation energies; this validates digitization scale but does not replace the reported values.",
            "An absolute complex epsilon_am*(omega,T) remains underdetermined until the high-frequency amorphous plateau is independently constrained. A one-Debye numerical bridge exists in code but requires that extra input explicitly.",
            "No v0.1.12 quantity is promoted automatically to TDGL eps_b or phase-local OAF parameters.",
        ],
    }

    outdir = ROOT / "outputs" / "v0112_same_state_bds"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"\nsaved: {path}")


if __name__ == "__main__":
    main()
