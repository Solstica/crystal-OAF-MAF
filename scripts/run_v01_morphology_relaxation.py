from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v01 import DIELECTRIC_CALIBRATION, PHASE_FRACTIONS
from pvdf_pf.calibration.dielectric import (
    effective_permittivity,
    infer_unknown_phase_permittivity,
    phase_eps_map,
)
from pvdf_pf.calibration.relaxation import (
    infer_debye_strength_from_real_epsilon,
    laminate_complex_permittivity,
    oaf_complex_permittivity,
)
from pvdf_pf.core.grid import Grid2D
from pvdf_pf.morphology.oriented import (
    periodic_winding_three_phase,
    winding_normal_angle_deg,
    wavy_winding_three_phase,
)
from pvdf_pf.morphology.three_phase import phase_fractions


def main() -> None:
    # Square cell used only for orientation/morphology sensitivity.  The physical
    # length scale is not calibrated in v0.1, so dx=dz=1 remains dimensionless.
    grid = Grid2D(nz=96, nx=96, dz=1.0, dx=1.0)

    crystal_eps = float(DIELECTRIC_CALIBRATION["crystal_eps_r"])
    maf_eps = float(DIELECTRIC_CALIBRATION["maf_eps_hypothesis"])
    target = float(DIELECTRIC_CALIBRATION["target_eps_eff"])

    # Transfer test: first infer the OAF effective permittivity from a perfectly
    # flat laminate with field parallel to lamellae, then hold that inferred value
    # fixed while rotating or wavifying the morphology.
    flat = periodic_winding_three_phase(
        grid,
        crystal_fraction=PHASE_FRACTIONS["crystal"],
        oaf_fraction=PHASE_FRACTIONS["oaf"],
        maf_fraction=PHASE_FRACTIONS["maf"],
        winding_z=1,
        winding_x=0,
    )
    oaf_eps_transfer = infer_unknown_phase_permittivity(
        flat,
        grid,
        target_eps_eff=target,
        known_phase_eps={"crystal": crystal_eps, "maf": maf_eps},
        unknown_phase="oaf",
        axis="x",  # x is parallel to flat lamellae when the normal is +z
        bracket=tuple(DIELECTRIC_CALIBRATION["oaf_fit_bracket"]),
    )
    phase_eps = {"crystal": crystal_eps, "oaf": oaf_eps_transfer, "maf": maf_eps}

    orientation_rows = []
    for winding_z, winding_x in [(1, 0), (2, 1), (1, 1), (1, 2), (0, 1)]:
        phase = periodic_winding_three_phase(
            grid,
            crystal_fraction=PHASE_FRACTIONS["crystal"],
            oaf_fraction=PHASE_FRACTIONS["oaf"],
            maf_fraction=PHASE_FRACTIONS["maf"],
            winding_z=winding_z,
            winding_x=winding_x,
        )
        eps_map = phase_eps_map(phase, phase_eps)
        eps_z, _ = effective_permittivity(eps_map, grid, axis="z")
        eps_x, _ = effective_permittivity(eps_map, grid, axis="x")
        orientation_rows.append(
            {
                "winding": [winding_z, winding_x],
                "lamellar_normal_angle_from_z_deg": winding_normal_angle_deg(grid, winding_z, winding_x),
                "phase_fractions_realized": phase_fractions(phase),
                "eps_eff_field_z": eps_z,
                "eps_eff_field_x": eps_x,
            }
        )

    waviness_rows = []
    for amplitude in [0.00, 0.03, 0.06, 0.10, 0.15]:
        phase = wavy_winding_three_phase(
            grid,
            crystal_fraction=PHASE_FRACTIONS["crystal"],
            oaf_fraction=PHASE_FRACTIONS["oaf"],
            maf_fraction=PHASE_FRACTIONS["maf"],
            winding_z=1,
            winding_x=0,
            waviness_amplitude=amplitude,
            waviness_mode_z=0,
            waviness_mode_x=2,
        )
        eps_map = phase_eps_map(phase, phase_eps)
        eps_z, _ = effective_permittivity(eps_map, grid, axis="z")
        eps_x, _ = effective_permittivity(eps_map, grid, axis="x")
        waviness_rows.append(
            {
                "waviness_amplitude_period": amplitude,
                "phase_fractions_realized": phase_fractions(phase),
                "eps_eff_field_z": eps_z,
                "eps_eff_field_x": eps_x,
            }
        )

    # Single-frequency identifiability scan.  tau is deliberately scanned rather
    # than fitted because epsilon'(10 Hz) alone cannot identify tau and Delta-eps.
    tau_rows = []
    for tau_s in [1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]:
        delta_eps = infer_debye_strength_from_real_epsilon(
            target_eps_real=target,
            frequency_hz=10.0,
            tau_s=tau_s,
            phase_fractions=PHASE_FRACTIONS,
            crystal_eps=crystal_eps,
            maf_eps=maf_eps,
            oaf_eps_fast=crystal_eps,
            orientation="parallel",
        )
        eps_oaf_10hz = complex(
            oaf_complex_permittivity(
                10.0,
                eps_fast=crystal_eps,
                delta_eps=delta_eps,
                tau_s=tau_s,
            )
        )
        eps_parallel, eps_normal = laminate_complex_permittivity(
            PHASE_FRACTIONS,
            {
                "crystal": complex(crystal_eps),
                "oaf": eps_oaf_10hz,
                "maf": complex(maf_eps),
            },
        )
        tau_rows.append(
            {
                "tau_s_hypothesis": tau_s,
                "delta_eps_oaf_required": delta_eps,
                "oaf_eps_real_10Hz": eps_oaf_10hz.real,
                "oaf_loss_epsilon_10Hz": -eps_oaf_10hz.imag,
                "laminate_parallel_eps_real_10Hz": eps_parallel.real,
                "laminate_parallel_loss_epsilon_10Hz": -eps_parallel.imag,
                "laminate_normal_eps_real_10Hz": eps_normal.real,
                "status": "INFERRED_CONDITIONAL",
            }
        )

    report = {
        "model_version": "v0.1-morphology-relaxation",
        "purpose": "test transferability of static phase permittivity and expose Debye parameter non-identifiability",
        "phase_fraction_target": PHASE_FRACTIONS,
        "static_transfer_hypothesis": {
            "crystal_eps_r": crystal_eps,
            "maf_eps_r": maf_eps,
            "oaf_eps_r_inferred_from_flat_parallel_case": oaf_eps_transfer,
            "target_eps_eff": target,
            "warning": "OAF value is an orientation-specific inverse parameter, not a measured phase constant.",
        },
        "orientation_scan": orientation_rows,
        "waviness_scan": waviness_rows,
        "debye_tau_scan": tau_rows,
        "interpretation_rules": [
            "Orientation/waviness results test whether one inverse phase-property set transfers across morphology.",
            "A single 10 Hz epsilon' value cannot determine both OAF Debye strength and relaxation time.",
            "The Debye scan is a conditional family; none of its tau values are promoted to physical parameters.",
            "Measured low-frequency epsilon is still kept separate from TDGL background eps_b to avoid double counting.",
        ],
    }

    outdir = ROOT / "outputs" / "v01_morphology_relaxation"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(json.dumps(report, indent=2))
    print(f"\nsaved: {path}")


if __name__ == "__main__":
    main()
