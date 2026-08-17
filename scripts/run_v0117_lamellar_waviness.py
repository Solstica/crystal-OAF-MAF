from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v0117 import (
    E_EXTERNAL_VPM,
    FROZEN_SOURCE,
    GRID,
    HUANG_SAXS_GEOMETRY_NM,
    OAF_PROFILE,
    PRONY,
    REPRESENTATIVE_STATE,
    RUI2022_DIR,
    STATIC_RELAXATION_FACTOR,
    STRUCTURE,
    V0117_RULES,
    WAVINESS_AMPLITUDES,
    WAVINESS_MODE_X,
    WAVINESS_MODE_Z,
)
from pvdf_pf.calibration.bds_full_spectrum import load_figure2_alpha_window, fit_relaxation
from pvdf_pf.calibration.prony import fit_positive_prony_from_cole_cole
from pvdf_pf.calibration.source_data import load_digitized_curve_csv
from pvdf_pf.core.grid import Grid2D
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


def waviness_normal_projection_metrics(
    grid: Grid2D,
    amplitude: float,
    *,
    mode_z: int,
    mode_x: int,
) -> dict[str, float]:
    """Analytic local interface-normal projection for the v0.1.17 phase coordinate."""
    z = np.arange(grid.nz, dtype=float)[:, None] / grid.nz
    x = np.arange(grid.nx, dtype=float)[None, :] / grid.nx
    phase = 2.0 * np.pi * (mode_z * z + mode_x * x)
    Lz = grid.nz * grid.dz
    Lx = grid.nx * grid.dx

    # u = x/Lx + A sin(2*pi*(mz*z/Lz + mx*x/Lx))
    du_dz = float(amplitude) * 2.0 * np.pi * mode_z * np.cos(phase) / Lz
    du_dx = 1.0 / Lx + float(amplitude) * 2.0 * np.pi * mode_x * np.cos(phase) / Lx
    denom = np.sqrt(du_dz**2 + du_dx**2)
    nz = np.divide(du_dz, denom, out=np.zeros_like(denom), where=denom > 0.0)
    return {
        "rms_abs_normal_projection_on_polarization": float(np.sqrt(np.mean(nz**2))),
        "mean_abs_normal_projection_on_polarization": float(np.mean(np.abs(nz))),
        "max_abs_normal_projection_on_polarization": float(np.max(np.abs(nz))),
    }


def through_origin_fit(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    denom = float(np.dot(x, x))
    if denom <= 0.0:
        return {"slope": 0.0, "r2": float("nan")}
    slope = float(np.dot(x, y) / denom)
    pred = slope * x
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot_origin = float(np.sum(y**2))
    r2 = float(1.0 - ss_res / ss_tot_origin) if ss_tot_origin > 0.0 else float("nan")
    return {"slope": slope, "r2_origin": r2}


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
    for amplitude in WAVINESS_AMPLITUDES:
        phase = wavy_winding_three_phase(
            grid,
            crystal_fraction=fc,
            oaf_fraction=fo,
            maf_fraction=fi,
            winding_z=0,
            winding_x=1,
            waviness_amplitude=float(amplitude),
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

        geometry = waviness_normal_projection_metrics(
            grid,
            float(amplitude),
            mode_z=WAVINESS_MODE_Z,
            mode_x=WAVINESS_MODE_X,
        )
        source_z_variation = float(
            np.sqrt(np.mean((P_frozen - np.roll(P_frozen, 1, axis=0)) ** 2))
        )
        rms_e_total = float(np.sqrt(np.mean(step.E_z**2 + step.E_x**2)))
        cases.append({
            "waviness_amplitude_fraction_period": float(amplitude),
            **geometry,
            "realized_phase_fractions": realized_phase_fractions(phase),
            "mean_P_frozen_beta_C_m2": phase_mean(P_frozen, phase, CRYSTAL),
            "mean_P_frozen_oaf_C_m2": phase_mean(P_frozen, phase, OAF),
            "mean_P_frozen_iaf_C_m2": phase_mean(P_frozen, phase, MAF),
            "rms_discrete_z_variation_P_frozen_C_m2": source_z_variation,
            "max_abs_Ez_Vpm": float(np.max(np.abs(step.E_z))),
            "max_abs_Ex_Vpm": float(np.max(np.abs(step.E_x))),
            "rms_Ez_Vpm": float(np.sqrt(np.mean(step.E_z**2))),
            "rms_Ex_Vpm": float(np.sqrt(np.mean(step.E_x**2))),
            "rms_E_total_Vpm": rms_e_total,
            "rms_E_total_over_Pscale_over_eps0": float(rms_e_total / field_scale),
            "mean_Ez_beta_Vpm": phase_mean(step.E_z, phase, CRYSTAL),
            "mean_Ez_oaf_Vpm": phase_mean(step.E_z, phase, OAF),
            "mean_Ez_iaf_Vpm": phase_mean(step.E_z, phase, MAF),
            "D_z_mean_face_C_m2": step.D_z_mean_face,
            "D_z_mean_cell_C_m2": step.D_z_mean_cell,
            "gauss_relative_residual": step.gauss_relative_residual,
            "constitutive_flux_mismatch": step.constitutive_flux_mismatch,
        })

    x = np.asarray(
        [c["rms_abs_normal_projection_on_polarization"] for c in cases], dtype=float
    )
    y_total = np.asarray([c["rms_E_total_Vpm"] for c in cases], dtype=float) / field_scale
    x2 = x**2
    y_z = np.asarray([c["rms_Ez_Vpm"] for c in cases], dtype=float) / field_scale

    total_fit = through_origin_fit(x, y_total)
    z_fit = through_origin_fit(x2, y_z)

    report = {
        "model_version": "v0.1.17",
        "stage": "Huang-aligned lamellar waviness sensitivity with frozen P_FE and generalized-Debye screening",
        "rules": V0117_RULES,
        "representative_bds_state": REPRESENTATIVE_STATE,
        "huang_saxs_geometry_nm": HUANG_SAXS_GEOMETRY_NM,
        "geometry_volume_fractions_from_thickness": {
            "beta_crystal": fc,
            "oaf": fo,
            "iaf": fi,
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
        "normalized_geometry_scaling": {
            "rms_total_field_vs_rms_abs_nz_through_origin": total_fit,
            "rms_Ez_vs_rms_nz_squared_through_origin": z_fit,
            "normalization_field_scale_Pmax_over_eps0_Vpm": field_scale,
        },
        "summary": {
            "flat_aligned_max_abs_Ez_Vpm": cases[0]["max_abs_Ez_Vpm"],
            "largest_waviness_amplitude": cases[-1]["waviness_amplitude_fraction_period"],
            "largest_waviness_rms_abs_nz": cases[-1][
                "rms_abs_normal_projection_on_polarization"
            ],
            "largest_waviness_rms_E_total_Vpm": cases[-1]["rms_E_total_Vpm"],
            "largest_waviness_max_abs_Ez_Vpm": cases[-1]["max_abs_Ez_Vpm"],
            "all_gauss_residuals_below_1e-7": bool(
                all(c["gauss_relative_residual"] < 1e-7 for c in cases)
            ),
            "tdgl_coupled": False,
            "ferroelectric_switching_enabled": False,
        },
        "interpretation": [
            "The zero-waviness Huang-aligned laminate is the v0.1.16 tangent-interface limit and should not create a frozen-source depolarization field.",
            "Sinusoidal waviness locally gives the beta/OAF/IAF interface normal a film-normal projection n_z. A discontinuity in film-normal P_FE then carries bound surface charge proportional to Delta(P_FE)*n_z.",
            "The sweep keeps the average lamellar orientation unchanged and keeps beta/OAF frozen-polarization amplitudes fixed. Field changes therefore isolate local geometric misorientation rather than a change in total frozen dipole moment.",
            "The waviness amplitudes are not fitted to Huang 2021. The available WAXD/SAXS text reports stronger MD orientation but does not provide a numerical azimuthal FWHM or a real-space interface-waviness distribution for this transfer.",
            "Absolute field values scale with the placeholder frozen-polarization amplitudes. The normalized geometry trends are the transferable result of this version.",
        ],
    }

    outdir = ROOT / "outputs" / "v0117_lamellar_waviness"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report["summary"], indent=2))
    print(json.dumps(report["normalized_geometry_scaling"], indent=2))
    print(f"saved: {path}")


if __name__ == "__main__":
    main()
