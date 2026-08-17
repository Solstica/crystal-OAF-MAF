"""Generate the 300 K dimensionless-rescaling audit for Ahluwalia et al. 2008."""
from __future__ import annotations

import json
from pathlib import Path

from pvdf_pf.literature.ahluwalia2008_dimensionless import (
    dimensionless_checkpoint_report,
)


OUT = Path("outputs/literature_reproduction/ahluwalia2008_dimensionless/report.json")


def main() -> None:
    report = {
        "benchmark": "Ahluwalia2008 Eq11-Eq12 300K dimensionless rescaling audit",
        "source": "Phys. Rev. B 78, 054110 (2008), DOI 10.1103/PhysRevB.78.054110",
        **dimensionless_checkpoint_report(),
        "reproducibility_guardrails": {
            "characteristic_length": (
                "The approximately 0.405 nm scale is inferred by inverting the source's "
                "reported 300 K z-noise rescaling relation. It is not a separately quoted "
                "source parameter and it is not the 2.16 nm TDGL grid spacing."
            ),
            "stochastic_trace": (
                "The article reports matched noise amplitudes and a 9 ps physical time "
                "mapping, but does not state an integration time step, random seed or a "
                "fully reproducible stochastic discretization. Exact point-by-point Fig. 5b "
                "trajectory reproduction is therefore not claimed."
            ),
            "next_solver_claim": (
                "A numerical TDGL implementation may be validated against coefficient "
                "identities, equilibrium statistics and reported physical scales, while its "
                "chosen time integrator/time step/seed must remain PROJECT_NUMERICAL_IMPLEMENTATION."
            ),
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
