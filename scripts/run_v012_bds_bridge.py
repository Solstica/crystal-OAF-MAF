from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v012 import BDS_PRIMARY_SOURCE, V012_RULES
from pvdf_pf.calibration.bds import debye_complex_permittivity, fit_debye_spectrum
from pvdf_pf.physics.dipolar import (
    complex_susceptibility_from_harmonic_history,
    discrete_debye_susceptibility,
    simulate_debye_history,
)


def main() -> None:
    # Purely synthetic self-test parameters. They are intentionally unrelated to
    # BOPVDF and must never be promoted to the physical parameter database.
    synthetic = {"eps_inf": 4.0, "delta_eps": 11.0, "tau_s": 2.0e-3}
    frequency = np.logspace(0.0, 5.0, 161)
    eps = debye_complex_permittivity(frequency, **synthetic)
    fit = fit_debye_spectrum(frequency, eps.real, -eps.imag)

    f_drive = 40.0
    steps_per_cycle = 400
    cycles = 40
    dt = 1.0 / (f_drive * steps_per_cycle)
    t = dt * np.arange(steps_per_cycle * cycles)
    E = 1.0e6 * np.cos(2.0 * np.pi * f_drive * t)
    P = simulate_debye_history(
        E,
        dt_s=dt,
        tau_s=synthetic["tau_s"],
        delta_eps=synthetic["delta_eps"],
    )
    chi_td = complex_susceptibility_from_harmonic_history(
        E,
        P,
        dt_s=dt,
        frequency_hz=f_drive,
        discard_fraction=0.5,
    )
    chi_zoh = discrete_debye_susceptibility(
        frequency_hz=f_drive,
        dt_s=dt,
        tau_s=synthetic["tau_s"],
        delta_eps=synthetic["delta_eps"],
    )
    chi_continuum = synthetic["delta_eps"] / (
        1.0 + 1j * 2.0 * np.pi * f_drive * synthetic["tau_s"]
    )

    report = {
        "model_version": "v0.1.2",
        "status": "NUMERICAL_BRIDGE_SELF_TEST",
        "primary_source_constraints": BDS_PRIMARY_SOURCE,
        "rules": V012_RULES,
        "synthetic_fit": {
            "input": synthetic,
            "recovered": fit.to_dict(),
            "warning": "Synthetic values are numerical verification only and are not BOPVDF parameters.",
        },
        "frequency_time_bridge": {
            "frequency_Hz": f_drive,
            "dt_s": dt,
            "continuum_Debye_delta_epsilon_complex": [
                float(chi_continuum.real),
                float(chi_continuum.imag),
            ],
            "zero_order_hold_delta_epsilon_complex": [float(chi_zoh.real), float(chi_zoh.imag)],
            "time_domain_delta_epsilon_complex": [float(chi_td.real), float(chi_td.imag)],
            "time_domain_vs_ZOH_absolute_error": float(abs(chi_td - chi_zoh)),
            "ZOH_vs_continuum_relative_error": float(
                abs(chi_zoh - chi_continuum) / abs(chi_continuum)
            ),
            "interpretation": (
                "The exact exponential update is exact for piecewise-constant (zero-order-hold) E. "
                "Its discrete harmonic transfer converges to continuum Debye response as dt -> 0."
            ),
        },
        "data_needed_for_physical_fit": [
            "BDS epsilon'(f,T) and epsilon''(f,T) for unpoled and poled BOPVDF",
            "paper-specific decomposition of amorphous/RAF/MOAF contributions",
            "numerical n(T), m_d(T), g(T), lambda(T), and/or mu_r(T) from source/SI where used",
            "validated RAF-to-OAF mapping",
            "physical TDGL time-scale calibration",
        ],
    }

    outdir = ROOT / "outputs" / "v012_bds_bridge"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"\nsaved: {path}")


if __name__ == "__main__":
    main()
