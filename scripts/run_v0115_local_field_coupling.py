from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v0115 import (
    DYNAMIC_REFERENCE_FREQUENCIES_HZ,
    E_EXTERNAL_VPM,
    GRID,
    PRONY,
    REPRESENTATIVE_STATE,
    RUI2022_DIR,
    STATIC_RELAXATION_FACTOR,
    STRUCTURE,
    V0115_RULES,
)
from pvdf_pf.calibration.bds_full_spectrum import load_figure2_alpha_window, fit_relaxation
from pvdf_pf.calibration.dielectric import effective_permittivity
from pvdf_pf.calibration.prony import fit_positive_prony_from_cole_cole
from pvdf_pf.calibration.source_data import load_digitized_curve_csv
from pvdf_pf.core.grid import Grid2D
from pvdf_pf.morphology.oriented import periodic_winding_three_phase
from pvdf_pf.morphology.three_phase import CRYSTAL, phase_fractions
from pvdf_pf.physics.coupled_local_field import (
    advance_self_consistent_generalized_debye,
    aligned_laminate_effective_permittivity,
    build_combined_amorphous_background_map,
)
from pvdf_pf.physics.dipolar import EPS0
from pvdf_pf.physics.generalized_debye import (
    GeneralizedDebyeBank,
    continuous_generalized_permittivity,
)


def exact_source_value(points, sample_state: str, temperature_C: float, quantity: str) -> float:
    values = [
        p.value for p in points
        if p.sample_state.lower() == sample_state.lower()
        and p.quantity == quantity
        and abs(p.temperature_C - temperature_C) < 1e-9
    ]
    if len(values) != 1:
        raise ValueError((sample_state, temperature_C, quantity, len(values)))
    return float(values[0])


def phase_mean(field: np.ndarray, phase: np.ndarray, mask_value: int | None) -> float:
    if mask_value is None:
        mask = phase != CRYSTAL
    else:
        mask = phase == mask_value
    return float(np.mean(field[mask]))


def main() -> None:
    sample = str(REPRESENTATIVE_STATE["sample_state"])
    T = float(REPRESENTATIVE_STATE["temperature_C"])

    spectra = load_figure2_alpha_window(RUI2022_DIR / "figure_2_bds_alpha_window_digitized.csv")
    static_points = load_digitized_curve_csv(RUI2022_DIR / "figure_3a_film_eps_digitized.csv")
    grouped = defaultdict(list)
    for p in spectra:
        grouped[(p.sample_state, p.temperature_C)].append(p)

    eps_film_static_source = exact_source_value(static_points, sample, T, "epsilon_c_film")
    cc = fit_relaxation(
        grouped[(sample, T)],
        epsilon_static_film=eps_film_static_source,
        model="cole-cole",
    )
    rep = fit_positive_prony_from_cole_cole(cc, **PRONY)
    bank = GeneralizedDebyeBank(
        epsilon_infinity=rep.epsilon_infinity_amorphous,
        tau_modes_s=rep.tau_modes_s,
        delta_epsilon_modes=rep.delta_epsilon_modes,
    )

    grid = Grid2D(**GRID)
    frac_kwargs = {
        "crystal_fraction": STRUCTURE["crystal_fraction"],
        "oaf_fraction": STRUCTURE["oaf_fraction"],
        "maf_fraction": STRUCTURE["maf_fraction"],
    }
    phase_parallel = periodic_winding_three_phase(
        grid, **frac_kwargs, winding_z=0, winding_x=1
    )
    phase_series = periodic_winding_three_phase(
        grid, **frac_kwargs, winding_z=1, winding_x=0
    )

    dt_static = float(STATIC_RELAXATION_FACTOR) * max(bank.tau_modes_s)
    P0 = np.zeros((bank.n_modes,) + grid.shape, dtype=float)

    geometry_results = {}
    for name, phase, is_parallel in (
        ("parallel_to_layers", phase_parallel, True),
        ("normal_to_layers", phase_series, False),
    ):
        eps_b, amorphous_mask = build_combined_amorphous_background_map(
            phase,
            bank,
            crystal_phase_id=CRYSTAL,
            epsilon_crystal_background=STRUCTURE["epsilon_crystal_background"],
        )
        step = advance_self_consistent_generalized_debye(
            P0,
            eps_b,
            amorphous_mask,
            grid,
            E_external_z=E_EXTERNAL_VPM,
            dt_s=dt_static,
            bank=bank,
        )

        eps_static_map = np.full(grid.shape, bank.epsilon_static_amorphous, dtype=float)
        eps_static_map[phase == CRYSTAL] = STRUCTURE["epsilon_crystal_background"]
        eps_eff_reference, _ = effective_permittivity(
            eps_static_map, grid, axis="z", E0=1.0
        )
        analytic = aligned_laminate_effective_permittivity(
            STRUCTURE["epsilon_crystal_background"],
            bank.epsilon_static_amorphous,
            crystal_fraction=STRUCTURE["crystal_fraction"],
            field_parallel_to_layers=is_parallel,
        )
        eps_eff_coupled = step.D_z_mean_face / (EPS0 * E_EXTERNAL_VPM)
        crystal_mean_E = phase_mean(step.E_z, phase, CRYSTAL)
        amorphous_mean_E = phase_mean(step.E_z, phase, None)

        geometry_results[name] = {
            "phase_fractions": phase_fractions(phase),
            "epsilon_eff_coupled_static_limit": eps_eff_coupled,
            "epsilon_eff_static_cell_reference": eps_eff_reference,
            "epsilon_eff_analytic_aligned_laminate": analytic.real,
            "relative_error_vs_cell_reference": (eps_eff_coupled - eps_eff_reference) / eps_eff_reference,
            "relative_error_vs_analytic": (eps_eff_coupled - analytic.real) / analytic.real,
            "mean_Ez_crystal_over_E0": crystal_mean_E / E_EXTERNAL_VPM,
            "mean_Ez_combined_amorphous_over_E0": amorphous_mean_E / E_EXTERNAL_VPM,
            "max_abs_Ez_over_E0": float(np.max(np.abs(step.E_z)) / E_EXTERNAL_VPM),
            "gauss_relative_residual": step.gauss_relative_residual,
            "constitutive_flux_mismatch": step.constitutive_flux_mismatch,
        }

    dynamic_reference = []
    for f in DYNAMIC_REFERENCE_FREQUENCIES_HZ:
        eps_am = complex(continuous_generalized_permittivity(bank, f))
        eps_parallel = aligned_laminate_effective_permittivity(
            STRUCTURE["epsilon_crystal_background"],
            eps_am,
            crystal_fraction=STRUCTURE["crystal_fraction"],
            field_parallel_to_layers=True,
        )
        eps_series = aligned_laminate_effective_permittivity(
            STRUCTURE["epsilon_crystal_background"],
            eps_am,
            crystal_fraction=STRUCTURE["crystal_fraction"],
            field_parallel_to_layers=False,
        )
        dynamic_reference.append({
            "frequency_Hz": float(f),
            "epsilon_amorphous_real": eps_am.real,
            "epsilon_amorphous_loss": -eps_am.imag,
            "epsilon_parallel_real": eps_parallel.real,
            "epsilon_parallel_loss": -eps_parallel.imag,
            "epsilon_series_real": eps_series.real,
            "epsilon_series_loss": -eps_series.imag,
            "series_crystal_field_amplitude_over_E0": abs(eps_series / STRUCTURE["epsilon_crystal_background"]),
            "series_amorphous_field_amplitude_over_E0": abs(eps_series / eps_am),
        })

    parallel_eps = geometry_results["parallel_to_layers"]["epsilon_eff_coupled_static_limit"]
    series_eps = geometry_results["normal_to_layers"]["epsilon_eff_coupled_static_limit"]
    report = {
        "model_version": "v0.1.15",
        "stage": "self-consistent generalized-Debye local-field coupling with frozen TDGL",
        "rules": V0115_RULES,
        "representative_state": REPRESENTATIVE_STATE,
        "source_and_constitutive": {
            "epsilon_film_static_source": eps_film_static_source,
            "epsilon_amorphous_static": bank.epsilon_infinity + bank.delta_epsilon_total,
            "epsilon_amorphous_infinity": bank.epsilon_infinity,
            "prony_active_modes": rep.active_mode_count,
            "prony_max_normalized_complex_error": rep.normalized_max_complex_error,
            "quasi_static_dt_s": dt_static,
        },
        "geometry_results": geometry_results,
        "dynamic_aligned_laminate_reference": dynamic_reference,
        "summary": {
            "parallel_geometry_recovers_source_static_film": abs(parallel_eps - eps_film_static_source) / eps_film_static_source < 5e-4,
            "parallel_epsilon_static": parallel_eps,
            "source_epsilon_static": eps_film_static_source,
            "series_epsilon_static": series_eps,
            "series_to_parallel_static_ratio": series_eps / parallel_eps,
            "series_crystal_field_amplification": geometry_results["normal_to_layers"]["mean_Ez_crystal_over_E0"],
            "series_amorphous_field_ratio": geometry_results["normal_to_layers"]["mean_Ez_combined_amorphous_over_E0"],
            "maximum_gauss_relative_residual": max(
                v["gauss_relative_residual"] for v in geometry_results.values()
            ),
            "tdgl_coupled": False,
            "constitutive_scope": "combined OAF+IAF amorphous response",
        },
        "interpretation": [
            "The generalized-Debye auxiliary polarization now enters Gauss' law self-consistently. One algebraically eliminated heterogeneous Poisson solve is sufficient per zero-order-hold time step.",
            "When the field lies parallel to the ideal lamellae, the self-consistent static limit reproduces the source two-phase film closure by construction and serves as a control geometry.",
            "When the field is normal to the same lamellae, dielectric contrast redistributes the local field and substantially changes the macroscopic response even though all phase constitutive parameters are unchanged.",
            "All non-crystal cells still carry the same combined-amorphous constitutive bank. OAF and IAF/MAF have not yet been separated by independent dynamic data.",
            "The ferroelectric TDGL order parameter remains frozen; v0.1.15 tests electrostatic feedback only and does not introduce an uncalibrated mapping between seconds and TDGL time.",
        ],
    }

    outdir = ROOT / "outputs" / "v0115_local_field_coupling"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report["summary"], indent=2))
    print("dynamic_reference:")
    print(json.dumps(dynamic_reference, indent=2))
    print(f"saved: {path}")


if __name__ == "__main__":
    main()
