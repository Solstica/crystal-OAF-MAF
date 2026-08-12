from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v014 import RUI2022_DIR
from configs.bopvdf_v016 import FIGURE5B_CSV, V016_RULES
from pvdf_pf.calibration.project_devitrification import build_regularized_oaf_state
from pvdf_pf.calibration.rui2022_state import build_rui2022_temperature_state


def main() -> None:
    temperatures = [-30, -20, -10, 0, 10, 20, 30, 40]
    rows = []
    for sample in ("unpoled", "poled"):
        for T in temperatures:
            frac = build_regularized_oaf_state(FIGURE5B_CSV, temperature_C=T)
            src = build_rui2022_temperature_state(
                RUI2022_DIR,
                temperature_C=T,
                sample_state=sample,
                mobility_boundary="l1",
            )
            rows.append(
                {
                    "sample_state": sample,
                    "temperature_C": T,
                    "q_devitrification": frac.q_devitrification,
                    "eta_ROAF": frac.roaf,
                    "eta_MOAF": frac.moaf,
                    "eta_IAF": frac.iaf,
                    "eta_crystal": frac.crystal,
                    "x_RAF_raw": frac.x_raf_raw,
                    "x_MAF_raw": frac.x_maf_raw,
                    "q_internal_disagreement": frac.q_internal_disagreement,
                    "epsilon_MOAF_source_target": src.epsilon_moaf,
                    "lambda_source": src.interaction_lambda,
                    "mu_r_l1_source_m2_per_Vs": src.rotational_mobility_m2_per_Vs,
                }
            )

    report = {
        "model_version": "v0.1.6",
        "stage": "temperature-driven project-regularized OAF devitrification",
        "rules": V016_RULES,
        "rows": rows,
        "guardrails": [
            "q(T) is a project regularization, not Rui-2022 Eq. S4.",
            "epsilon_MOAF is not copied into TDGL eps_b.",
            "mu_r is not converted into tau without an independently validated relation.",
            "spatial ROAF/MOAF placement still requires an external mobility_score field.",
        ],
        "status": "REGULARIZED_OAF_TEMPERATURE_STATE_READY",
    }

    outdir = ROOT / "outputs" / "v016_devitrification_state"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"\nsaved: {path}")


if __name__ == "__main__":
    main()
