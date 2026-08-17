from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v0118 import (
    E_EXTERNAL_VPM,
    FROZEN_SOURCE,
    GRID,
    HUANG_SAXS_GEOMETRY_NM,
    OAF_PROFILE,
    ORIENTATION_PROVENANCE,
    PRONY,
    REPRESENTATIVE_STATE,
    RMS_NZ_SCAN,
    RUI2022_DIR,
    SENSITIVITY_BOUND_RMS_NZ,
    SOURCE_CONSTRAINED_RMS_NZ,
    STATIC_RELAXATION_FACTOR,
    STRUCTURE,
    V0118_RULES,
    WAVINESS_MODE_X,
    WAVINESS_MODE_Z,
)
from pvdf_pf.calibration.bds_full_spectrum import load_figure2_alpha_window, fit_relaxation
from pvdf_pf.calibration.prony import fit_positive_prony_from_cole_cole
from pvdf_pf.calibration.source_data import load_digitized_curve_csv
from pvdf_pf.core.grid import Grid2D
from pvdf_pf.morphology.orientation_constraint import (
    amplitude_for_target_rms_nz,
    sinusoidal_normal_projection_metrics,
)
from pvdf_pf.morphology.oriented import wavy_winding_three_phase
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


def through_origin_fit(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    denom = float(np.dot(x, x))
    if denom <= 0.0:
        return {"slope": 0.0, "r2_origin": float("nan")}
    slope = float(np.dot(x, y) / denom)
    pred = slope * x
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot_origin = float(np.sum(y**2))
    r2 = float(1.0 - ss_res / ss_tot_origin) if ss_tot_origin > 0.0 else float("nan")
    return {"slope": slope, "r2_origin": r2}


def provenance_label(target: float) -> str:
    if target == 0.0:
        return "IDEAL_TANGENT_CONTROL"
    if any(np.isclose(target, x) for x in SOURCE_CONSTRAINED_RMS_NZ):
        return "PROJECT_IMAGE_DERIVED_CONSTRAINT"
    if any(np.isclose(target, x) for x in SENSITIVITY_BOUND_RMS_NZ):
        return "PLACEHOLDER_SENSITIVITY_BOUND"
    raise ValueError(f"unclassified target RMS n_z: {target}")


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
    for target_rms_nz in RMS_NZ_SCAN:
        amplitude = amplitude_for_target_rms_nz(
            grid,
            float(target_rms_nz),
            mode_z=WAVINESS_MODE_Z,
            mode_x=WAVINESS_MODE_X,
        )
        geometry = sinusoidal_normal_projection_metrics(
            grid,
            amplitude,
            mode_z=WAVINESS_MODE_Z,
            mode_x=WAVINESS_MODE_X,
        )

        phase = wavy_winding_three_phase(
            grid,
            crystal_fraction=fc,
            oaf_fraction=fo,
            maf_fraction=fi,
            winding_z=0,
            winding_x=1,
            waviness_amplitude=amplitude,
            waviness_mode_z=WAVINESS_MODE_Z,
            waviness_mode_x=WAVINESS_MODE_X,
        )
        eps_b, amorphous_mask = build_combined_amorphous_background_map(
            phase,
            bank,
            crystal_phase_id=CRYSTAL,
            epsilon_crystal_background=STRUCTURE["epsilon_crystal_background"],
        )
        P_frozen = build_lamellar_frozen_polarization_z(
            phase,
            grid,
            crystal_fraction=fc,
            oaf_fraction=fo,
            p_beta_C_m2=FROZEN_SOURCE["p_beta_C_m2"],
            p_oaf_mean_C_m2=FROZEN_SOURCE["p_oaf_mean_C_m2"],
            winding_z=0,
            winding_x=1,
            oaf_profile=OAF_PROFILE,
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

        rms_e_total = float(np.sqrt(np.mean(step.E_z**2 + step.E_x**2)))
        cases.append({
            "target_rms_n_ND": float(target_rms_nz),
            "realized_rms_n_ND": geometry["rms_nz"],
            "orientation_provenance": provenance_label(float(target_rms_nz)),
            "equivalent_internal_waviness_amplitude_fraction_period": amplitude,
            "mean_abs_n_ND": geometry["mean_abs_nz"],
            "max_abs_n_ND": geometry["max_abs_nz"],
            "realized_phase_fractions": realized_phase_fractions(phase),
            "mean_P_frozen_beta_C_m2": phase_mean(P_frozen, phase, CRYSTAL),
            "mean_P_frozen_oaf_C_m2": phase_mean(P_frozen, phase, OAF),
            "mean_P_frozen_iaf_C_m2": phase_mean(P_frozen, phase, MAF),
            "rms_Ez_Vpm": float(np.sqrt(np.mean(step.E_z**2))),
            "rms_Ex_Vpm": float(np.sqrt(np.mean(step.E_x**2))),
            "rms_E_total_Vpm": rms_e_total,
            "rms_E_total_over_Pscale_over_eps0": float(rms_e_total / field_scale),
            "max_abs_Ez_Vpm": float(np.max(np.abs(step.E_z))),
            "max_abs_Ex_Vpm": float(np.max(np.abs(step.E_x))),
            "D_z_mean_face_C_m2": step.D_z_mean_face,
            "gauss_relative_residual": step.gauss_relative_residual,
            "constitutive_flux_mismatch": step.constitutive_flux_mismatch,
        })

    constrained = [c for c in cases if c["orientation_provenance"] == "PROJECT_IMAGE_DERIVED_CONSTRAINT"]
    x = np.asarray([c["realized_rms_n_ND"] for c in constrained])
    y = np.asarray([c["rms_E_total_over_Pscale_over_eps0"] for c in constrained])
    constrained_fit = through_origin_fit(x, y)

    report = {
        "model_version": "v0.1.18",
        "stage": "Huang-SAXS image-constrained lamellar-normal spread scan",
        "rules": V0118_RULES,
        "orientation_provenance": ORIENTATION_PROVENANCE,
        "representative_bds_state": REPRESENTATIVE_STATE,
        "huang_saxs_geometry_nm": HUANG_SAXS_GEOMETRY_NM,
        "frozen_source_sensitivity": FROZEN_SOURCE,
        "constitutive": {
            "epsilon_film_static_source": eps_film_static_source,
            "epsilon_combined_amorphous_static": bank.epsilon_static,
            "epsilon_combined_amorphous_infinity": bank.epsilon_infinity,
            "epsilon_beta_background": STRUCTURE["epsilon_crystal_background"],
            "prony_active_modes": rep.active_mode_count,
            "quasi_static_dt_s": dt_static,
        },
        "cases": cases,
        "source_constrained_geometry_fit": {
            "normalized_rms_total_field_vs_rms_n_ND": constrained_fit,
            "normalization_field_scale_Pmax_over_eps0_Vpm": field_scale,
        },
        "summary": {
            "number_source_constrained_cases": len(constrained),
            "source_constrained_rms_n_ND_min": min(c["realized_rms_n_ND"] for c in constrained),
            "source_constrained_rms_n_ND_max": max(c["realized_rms_n_ND"] for c in constrained),
            "source_constrained_normalized_rms_field_min": min(c["rms_E_total_over_Pscale_over_eps0"] for c in constrained),
            "source_constrained_normalized_rms_field_max": max(c["rms_E_total_over_Pscale_over_eps0"] for c in constrained),
            "maximum_target_realization_error": max(abs(c["target_rms_n_ND"] - c["realized_rms_n_ND"]) for c in cases),
            "all_gauss_residuals_below_1e-7": bool(all(c["gauss_relative_residual"] < 1e-7 for c in cases)),
            "tdgl_coupled": False,
            "ferroelectric_switching_enabled": False,
        },
        "interpretation": [
            "v0.1.18 reports the geometry by RMS film-normal interface-normal projection rather than by the internal sinusoidal amplitude.",
            "The interval 0.05-0.20 is deliberately wider than the raster-image central estimate and is tagged PROJECT_IMAGE_DERIVED_CONSTRAINT; Huang et al. did not directly report these RMS values.",
            "The 0.25 and 0.30 cases are sensitivity bounds and must not be presented as source-constrained BOPVDF texture.",
            "Absolute field magnitudes remain proportional to the normalized frozen beta/OAF polarization amplitudes. The current quantitative output is therefore the normalized field response to a constrained geometry descriptor.",
        ],
    }

    outdir = ROOT / "outputs" / "v0118_orientation_constrained_scan"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report["summary"], indent=2))
    print(json.dumps(report["source_constrained_geometry_fit"], indent=2))
    print(f"saved: {path}")


if __name__ == "__main__":
    main()
