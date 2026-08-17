"""Generate the Ahluwalia 2008 MD -> LGD/TDGL parameter-transfer audit."""
from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from pvdf_pf.literature.ahluwalia2008 import (
    GRADIENT_CHECKPOINT,
    KINETIC_CHECKPOINT,
    TABLE_I_PRINTED,
    TABLE_II_PRINTED,
    derive_lgd_from_table_i,
    kinetic_ratios,
    table_ii_rounding_audit,
)


OUT = Path("outputs/literature_reproduction/ahluwalia2008_multiscale/report.json")


def main() -> None:
    report = {
        "benchmark": "Ahluwalia2008 MD-to-LGD/TDGL multiscale parameter-transfer audit",
        "source": "Phys. Rev. B 78, 054110 (2008), DOI 10.1103/PhysRevB.78.054110",
        "direct_source": {
            "table_I_MD_observables": asdict(TABLE_I_PRINTED),
            "table_II_LGD_parameters": asdict(TABLE_II_PRINTED),
            "gradient_checkpoint": asdict(GRADIENT_CHECKPOINT),
            "kinetic_checkpoint": asdict(KINETIC_CHECKPOINT),
        },
        "project_rederived_from_rounded_source_values": {
            "LGD_from_Eqs_2_3_and_printed_Table_I": asdict(derive_lgd_from_table_i()),
            "rounding_audit_against_printed_Table_II": table_ii_rounding_audit(),
            "kinetic_ratios_from_Eq_11_noise_amplitudes": kinetic_ratios(),
        },
        "provenance_guardrails": {
            "table_II_authority": (
                "The printed Table II remains authoritative. Re-deriving it from the rounded "
                "Table I necessarily produces percent-level differences and must not replace "
                "the directly printed coefficients."
            ),
            "K3": (
                "K1 and K2 are obtained from MD-estimated 180-degree domain-wall widths; "
                "K3 is not MD-calibrated and is set equal to K1=K2 by the source authors for "
                "computational convenience."
            ),
            "scope": (
                "The source intentionally excludes static defects and assumes perfect "
                "compensation when studying intrinsic switching. These are source-model "
                "assumptions, not generic properties of semicrystalline PVDF."
            ),
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
