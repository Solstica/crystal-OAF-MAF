from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v0116 import (
    E_EXTERNAL_VPM,
    FROZEN_SOURCE,
    GRID,
    HUANG_SAXS_GEOMETRY_NM,
    OAF_PROFILES,
    ORIENTATIONS,
    PRONY,
    REPRESENTATIVE_STATE,
    RUI2022_DIR,
    STATIC_RELAXATION_FACTOR,
    STRUCTURE,
    V0116_RULES,
)
from pvdf_pf.calibration.bds_full_spectrum import load_figure2_alpha_window, fit_relaxation
from pvdf_pf.calibration.prony import fit_positive_prony_from_cole_cole
from pvdf_pf.calibration.source_data import load_digitized_curve_csv
from pvdf_pf.core.grid import Grid2D
from pvdf_pf.morphology.oriented import periodic_winding_three_phase, winding_normal_angle_deg
from pvdf_pf.morphology.three_phase import CRYSTAL, OAF, MAF
from pvdf_pf.physics.coupled_local_field import build_combined_amorphous_background_map
from pvdf_pf.physics.dipolar import EPS0
from pvdf_pf.physics.frozen_source import (
    advance_generalized_debye_with_frozen_source,
    build_lamellar_frozen_polarization_z,
)
from pvdf_pf.physics.generalized_debye import GeneralizedDebyeBank


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


def phase_mean(field: np.ndarray, phase: np.ndarray, phase_id: int) -> float:
    return float(np.mean(np.asarray(field)[phase == phase_id]))


def realized_phase_fractions(phase: np.ndarray) -> dict[str, float]:
    n = float(phase.size)
    return {
        "beta_crystal": float(np.count_nonzero(phase == CRYSTAL) / n),
        "oaf": float(np.count_nonzero(phase == OAF) / n),
        "iaf": float(np.count_nonzero(phase == MAF) / n),
    }


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
    fc = float(STRUCTURE["crystal_fraction_volume"])
    fo = float(STRUCTURE["oaf_fraction_volume"])
    fi = float(STRUCTURE["iaf_fraction_volume"])
    dt_static = float(STATIC_RELAXATION_FACTOR) * max(bank.tau_modes_s)
    P0 = np.zeros((bank.n_modes,) + grid.shape, dtype=float)
    field_scale = max(
        abs(float(FROZEN_SOURCE["p_beta_C_m2"])),
        abs(float(FROZEN_SOURCE["p_oaf_mean_C_m2"])),
    ) / EPS0

    cases = []
    for orientation in ORIENTATIONS:
        wz = int(orientation["winding_z"])
        wx = int(orientation["winding_x"])
        phase = periodic_winding_three_phase(
            grid,
            crystal_fraction=fc,
            oaf_fraction=fo,
            maf_fraction=fi,
            winding_z=wz,
            winding_x=wx,
        )
        eps_b, amorphous_mask = build_combined_amorphous_background_map(
            phase,
            bank,
            crystal_phase_id=CRYSTAL,
            epsilon_crystal_background=STRUCTURE["epsilon_crystal_background"],
        )
        angle = winding_normal_angle_deg(grid, wz, wx)
        nz_projection = abs(float(np.cos(np.deg2rad(angle))))

        for profile in OAF_PROFILES:
            P_frozen = build_lamellar_frozen_polarization_z(
                phase,
                grid,
                crystal_fraction=fc,
                oaf_fraction=fo,
                p_beta_C_m2=FROZEN_SOURCE["p_beta_C_m2"],
                p_oaf_mean_C_m2=FROZEN_SOURCE["p_oaf_mean_C_m2"],
                winding_z=wz,
                winding_x=wx,
                oaf_profile=profile,
                decay_length_fraction_of_half_oaf=FROZEN_SOURCE[
                    "oaf_decay_length_fraction_of_half_oaf"
                ],
            )
            step = advance_generalized_debye_with_frozen_source(
                P0,
                eps_b,
                amorphous_mask,
                P_frozen,
                grid,
                E_external_z=E_EXTERNAL_VPM,
                dt_s=dt_static,
                bank=bank,
            )

            source_z_variation = float(
                np.sqrt(np.mean((P_frozen - np.roll(P_frozen, 1, axis=0)) ** 2))
            )
            cases.append({
                "orientation": orientation["name"],
                "winding_z": wz,
                "winding_x": wx,
                "lamellar_normal_angle_from_polarization_deg": angle,
                "abs_normal_projection_on_polarization": nz_projection,
                "oaf_profile": profile,
                "realized_phase_fractions": realized_phase_fractions(phase),
                "mean_P_frozen_beta_C_m2": phase_mean(P_frozen, phase, CRYSTAL),
                "mean_P_frozen_oaf_C_m2": phase_mean(P_frozen, phase, OAF),
                "mean_P_frozen_iaf_C_m2": phase_mean(P_frozen, phase, MAF),
                "mean_P_frozen_cell_C_m2": float(np.mean(P_frozen)),
                "rms_discrete_z_variation_P_frozen_C_m2": source_z_variation,
                "max_abs_Ez_Vpm": float(np.max(np.abs(step.E_z))),
                "rms_Ez_Vpm": float(np.sqrt(np.mean(step.E_z**2))),
                "mean_Ez_Vpm": float(np.mean(step.E_z)),
                "max_abs_Ez_over_Pscale_over_eps0": float(np.max(np.abs(step.E_z)) / field_scale),
                "mean_Ez_beta_Vpm": phase_mean(step.E_z, phase, CRYSTAL),
                "mean_Ez_oaf_Vpm": phase_mean(step.E_z, phase, OAF),
                "mean_Ez_iaf_Vpm": phase_mean(step.E_z, phase, MAF),
                "D_z_mean_face_C_m2": step.D_z_mean_face,
                "D_z_mean_cell_C_m2": step.D_z_mean_cell,
                "gauss_relative_residual": step.gauss_relative_residual,
                "constitutive_flux_mismatch": step.constitutive_flux_mismatch,
            })

    by_key = {(c["orientation"], c["oaf_profile"]): c for c in cases}
    tangent = by_key[("huang_aligned_tangent", "uniform")]
    normal = by_key[("normal_to_layers_control", "uniform")]
    decay_normal = by_key[("normal_to_layers_control", "interface_decay")]

    report = {
        "model_version": "v0.1.16",
        "stage": "frozen ferroelectric source plus generalized-Debye local-field orientation audit",
        "rules": V0116_RULES,
        "representative_bds_state": REPRESENTATIVE_STATE,
        "huang_saxs_geometry_nm": HUANG_SAXS_GEOMETRY_NM,
        "geometry_volume_fractions_from_thickness": {
            "beta_crystal": fc,
            "oaf": fo,
            "iaf": fi,
        },
        "source_fraction_metadata": {
            "beta_weight_fraction_waxd": STRUCTURE["source_beta_weight_fraction_waxd"],
            "minimum_oaf_weight_fraction_source_model": STRUCTURE[
                "source_min_oaf_weight_fraction_model"
            ],
        },
        "frozen_source_sensitivity": FROZEN_SOURCE,
        "constitutive": {
            "epsilon_film_static_source": eps_film_static_source,
            "epsilon_combined_amorphous_static": bank.epsilon_static,
            "epsilon_combined_amorphous_infinity": bank.epsilon_infinity,
            "epsilon_beta_background": STRUCTURE["epsilon_crystal_background"],
            "prony_active_modes": rep.active_mode_count,
            "prony_max_normalized_complex_error": rep.normalized_max_complex_error,
            "quasi_static_dt_s": dt_static,
        },
        "cases": cases,
        "summary": {
            "v0115_regression_fixed_before_v0116": True,
            "huang_aligned_uniform_max_abs_Ez_Vpm": tangent["max_abs_Ez_Vpm"],
            "huang_aligned_uniform_source_z_variation_C_m2": tangent[
                "rms_discrete_z_variation_P_frozen_C_m2"
            ],
            "normal_control_uniform_max_abs_Ez_Vpm": normal["max_abs_Ez_Vpm"],
            "normal_control_decay_max_abs_Ez_Vpm": decay_normal["max_abs_Ez_Vpm"],
            "normal_decay_to_uniform_peak_field_ratio": decay_normal["max_abs_Ez_Vpm"]
            / normal["max_abs_Ez_Vpm"],
            "ideal_huang_lamellae_frozen_source_generates_depolarization_field": tangent[
                "max_abs_Ez_Vpm"
            ] > 1e-6 * field_scale,
            "tdgl_coupled": False,
            "ferroelectric_switching_enabled": False,
        },
        "interpretation": [
            "The Huang SAXS geometry is transferred by thickness ratios: 5.78 nm beta crystal, 1.51 nm OAF on each side, and 3.00 nm IAF in an 11.8 nm period. The WAXD beta fraction 0.52 is kept as a weight-fraction datum rather than voxel occupancy.",
            "For the experimentally motivated aligned orientation, the lamellar normal is along MD while the poled polarization is along the film normal. The prescribed scalar Pz is tangent to ideal flat crystal/OAF interfaces, so its discontinuity occurs along x and does not create div(Pz ez).",
            "Tilted and normal-to-layer controls intentionally give the morphology normal a z component. They therefore create bound-charge source terms and local depolarization fields from the same frozen polarization amplitudes.",
            "Uniform and interface-decaying OAF profiles have the same discrete OAF mean. Differences between them isolate spatial allocation rather than changing the total OAF dipole moment.",
            "The frozen beta/OAF amplitudes are normalized sensitivity values. Available bulk remanent polarization does not uniquely identify phase-resolved remanent polarization, so v0.1.16 does not claim an absolute local P_FE calibration.",
        ],
    }

    outdir = ROOT / "outputs" / "v0116_frozen_source_orientation"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report["summary"], indent=2))
    print(f"saved: {path}")


if __name__ == "__main__":
    main()
