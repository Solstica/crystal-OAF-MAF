"""Generate source-faithful algebraic checkpoints for Su et al. 2022 PVDF Eq. (4)."""
from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from pvdf_pf.literature.su2022 import SU2022_PVDF, z_axis_global_minimum


OUT = Path("outputs/literature_reproduction/su2022_pvdf_landau/report.json")
CHECK_TEMPERATURES_K = (300.0, 315.0, 330.0)


def main() -> None:
    report = {
        "benchmark": "Su2022 beta-PVDF homogeneous Landau Eq. (4)",
        "source": "Nature Communications 13, 4867 (2022), DOI 10.1038/s41467-022-32518-3",
        "direct_source": {
            "pvdf_supplementary_table_3": asdict(SU2022_PVDF),
            "bulk_energy": (
                "f = alpha1 Px^2 + alpha2 Py^2 + alpha3(T) Pz^2 "
                "+ alpha33 Pz^4 + alpha333 Pz^6"
            ),
            "electrostatic_equilibrium": "div(eps0*eps_b*E + P) = 0",
            "reported_domain_simulation": {
                "domain_nm": [512.0, 512.0, 512.0],
                "grid": [128, 128, 128],
                "spacing_nm": 4.0,
                "applied_field_V_m": 1.2e5,
                "periodic_fields": True,
                "software": "MuPRO Ferroelectric module",
            },
        },
        "project_algebra_checkpoints_at_declared_temperatures": {
            str(int(t)): {
                "temperature_K": t,
                "alpha3_J_m_C2": SU2022_PVDF.alpha3(t),
                "nonnegative_z_axis_global_minimum": z_axis_global_minimum(t),
            }
            for t in CHECK_TEMPERATURES_K
        },
        "reproducibility_guardrails": {
            "temperature": (
                "The provided source files print the temperature-dependent coefficient law "
                "but do not identify a unique temperature for all phase-field runs. The "
                "300/315/330 K values above are declared algebra checks, not source results."
            ),
            "eps_b": (
                "The paper varies background dielectric constant with MXene loading. No "
                "project default is assigned in the homogeneous benchmark."
            ),
            "full_domain_solver": (
                "The total phase-field energy includes a gradient term, but the provided "
                "article/SI do not tabulate a PVDF gradient coefficient or kinetic coefficient. "
                "Therefore this report reproduces Eq. (4)/Table 3 algebra and Eq. (5) structure, "
                "not the full published domain morphology."
            ),
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
