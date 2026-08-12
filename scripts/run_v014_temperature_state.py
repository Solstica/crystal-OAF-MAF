from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v014 import MOBILITY_BOUNDARIES, RUI2022_DIR, SAMPLE_STATES, V014_RULES
from pvdf_pf.calibration.rui2022_state import build_rui2022_temperature_state


def main() -> None:
    temperatures = [-30, -20, -10, 0, 10, 20, 30, 40]
    states = []
    for sample in SAMPLE_STATES:
        for boundary in MOBILITY_BOUNDARIES:
            for T in temperatures:
                states.append(
                    build_rui2022_temperature_state(
                        RUI2022_DIR,
                        temperature_C=T,
                        sample_state=sample,
                        mobility_boundary=boundary,
                    ).to_dict()
                )

    report = {
        "model_version": "v0.1.4",
        "stage": "temperature-dependent source constitutive state",
        "rules": V014_RULES,
        "n_states": len(states),
        "states": states,
        "readiness": {
            "temperature_dependent_source_state": True,
            "epsilon_MOAF_target": True,
            "ROAF_MOAF_volume_fraction_field": False,
            "reason_fraction_field_blocked": (
                "Uploaded SI states that xRAF(T)/xMAF(T) are taken from main-text Figure 5B; "
                "those curves are not present in the SI digitization."
            ),
            "mu_r_to_tau": False,
            "TDGL_seconds_mapping": False,
        },
    }

    outdir = ROOT / "outputs" / "v014_temperature_state"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"\nsaved: {path}")


if __name__ == "__main__":
    main()
