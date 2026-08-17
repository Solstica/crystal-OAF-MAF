from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v0119 import (
    E_EXTERNAL_VPM,
    FROZEN_SOURCE,
    GRID,
    HUANG_SAXS_GEOMETRY_NM,
    OAF_PROFILE,
    ORIENTATION_SPECTRA,
    PRONY,
    REPRESENTATIVE_STATE,
    RUI2022_DIR,
    STATIC_RELAXATION_FACTOR,
    STRUCTURE,
    TARGET_RMS_N_ND,
    V0119_RULES,
)
from pvdf_pf.calibration.bds_full_spectrum import load_figure2_alpha_window, fit_relaxation
from pvdf_pf.calibration.prony import fit_positive_prony_from_cole_cole
from pvdf_pf.calibration.source_data import load_digitized_curve_csv
from pvdf_pf.core.grid import Grid2D
from pvdf_pf.morphology.orientation_ensemble import (
    HarmonicMode,
    scale_for_target_rms_normal,
    spectral_orientation_metrics,
    spectral_three_phase,
)
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


def parse_modes(spec: dict) -> tuple[HarmonicMode, ...]:
    return tuple(
        HarmonicMode(mode_z=int(k), coefficient=float(c), phase_rad=float(phi))
        for k, c, phi in spec["modes"]
    )


def phase_fraction_dict(phase: np.ndarray) -> dict[str, float]:
    n = float(phase.size)
    return {
        "beta": float(np.count_nonzero(phase == CRYSTAL) / n),
        "oaf": float(np.count_nonzero(phase == OAF) / n),
        "iaf": float(np.count_nonzero(phase == MAF) / n),
    }


def phase_rms(field: np.ndarray, phase: np.ndarray, phase_id: int) -> float:
    values = np.asarray(field, dtype=float)[phase == phase_id]
    return float(np.sqrt(np.mean(values**2)))


def normalized_source_divergence_rms(Pz: np.ndarray, grid, p_scale: float) -> float:
    """Discrete Pz-divergence diagnostic in inverse-grid-length units.

    The quantity is normalized by the frozen-polarization scale.  It is an internal
    morphology diagnostic, not a physical charge density because the current grid
    spacing is a normalized electrostatic audit length.
    """
    p = np.asarray(Pz, dtype=float) / float(p_scale)
    p_face_p = 0.5 * (p + np.roll(p, -1, axis=0))
    p_face_m = 0.5 * (p + np.roll(p, 1, axis=0))
    div = (p_face_p - p_face_m) / float(grid.dz)
    return float(np.sqrt(np.mean(div**2)))


def relative_spread(values: list[float]) -> dict[str, float]:
    x = np.asarray(values, dtype=float)
    mean = float(np.mean(x))
    return {
        "min": float(np.min(x)),
        "max": float(np.max(x)),
        "mean": mean,
        "range_over_mean": float((np.max(x) - np.min(x)) / mean) if mean != 0.0 else float("nan"),
        "coefficient_of_variation": float(np.std(x) / mean) if mean != 0.0 else float("nan"),
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
    p_scale = max(
        abs(float(FROZEN_SOURCE["p_beta_C_m2"])),
        abs(float(FROZEN_SOURCE["p_oaf_mean_C_m2"])),
    )
    e_scale = p_scale / EPS0

    cases = []
    for target in TARGET_RMS_N_ND:
        for spec in ORIENTATION_SPECTRA:
            modes = parse_modes(spec)
            scale = scale_for_target_rms_normal(grid, modes, float(target))
            metrics = spectral_orientation_metrics(grid, modes, scale)
            phase, n_nd = spectral_three_phase(
                grid,
                crystal_fraction=fc,
                oaf_fraction=fo,
                iaf_fraction=fi,
                modes=modes,
                scale=scale,
            )
            eps_b, amorphous_mask = build_combined_amorphous_background_map(
                phase,
                bank,
                crystal_phase_id=CRYSTAL,
                epsilon_crystal_background=STRUCTURE["epsilon_crystal_background"],
            )
            # Uniform OAF source is independent of the internal laminate coordinate,
            # so the existing phase-label builder remains valid for this ensemble.
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

            Ez_corr = step.E_z - float(E_EXTERNAL_VPM)
            E_mag_corr = np.sqrt(Ez_corr**2 + step.E_x**2)
            abs_e = np.abs(E_mag_corr) / e_scale
            cases.append({
                "target_rms_n_ND": float(target),
                "spectrum_name": spec["name"],
                "spectrum_status": spec["status"],
                "modes": [
                    {"mode_z": m.mode_z, "coefficient": m.coefficient, "phase_rad": m.phase_rad}
                    for m in modes
                ],
                "spectrum_scale": scale,
                "orientation_metrics": metrics,
                "direct_rms_check_from_n_field": float(np.sqrt(np.mean(n_nd**2))),
                "phase_fractions": phase_fraction_dict(phase),
                "normalized_discrete_source_divergence_rms": normalized_source_divergence_rms(P_frozen, grid, p_scale),
                "normalized_corrector_field": {
                    "rms": float(np.sqrt(np.mean(E_mag_corr**2)) / e_scale),
                    "mean": float(np.mean(E_mag_corr) / e_scale),
                    "p90": float(np.quantile(abs_e, 0.90)),
                    "p95": float(np.quantile(abs_e, 0.95)),
                    "p99": float(np.quantile(abs_e, 0.99)),
                    "max": float(np.max(abs_e)),
                    "rms_beta": phase_rms(E_mag_corr / e_scale, phase, CRYSTAL),
                    "rms_oaf": phase_rms(E_mag_corr / e_scale, phase, OAF),
                    "rms_iaf": phase_rms(E_mag_corr / e_scale, phase, MAF),
                },
                "gauss_relative_residual": step.gauss_relative_residual,
                "constitutive_flux_mismatch": step.constitutive_flux_mismatch,
            })

    comparisons = []
    for target in TARGET_RMS_N_ND:
        subset = [c for c in cases if np.isclose(c["target_rms_n_ND"], target)]
        comparisons.append({
            "target_rms_n_ND": float(target),
            "n_spectra": len(subset),
            "realized_rms_n_ND_spread": relative_spread([
                c["orientation_metrics"]["rms_n_ND"] for c in subset
            ]),
            "normalized_field_rms_spread": relative_spread([
                c["normalized_corrector_field"]["rms"] for c in subset
            ]),
            "normalized_field_p99_spread": relative_spread([
                c["normalized_corrector_field"]["p99"] for c in subset
            ]),
            "normalized_field_max_spread": relative_spread([
                c["normalized_corrector_field"]["max"] for c in subset
            ]),
            "source_divergence_rms_spread": relative_spread([
                c["normalized_discrete_source_divergence_rms"] for c in subset
            ]),
        })

    report = {
        "model_version": "v0.1.19",
        "stage": "same-RMS lamellar-normal orientation-spectrum sensitivity",
        "rules": V0119_RULES,
        "representative_bds_state": REPRESENTATIVE_STATE,
        "huang_saxs_geometry_nm": HUANG_SAXS_GEOMETRY_NM,
        "frozen_source_sensitivity": FROZEN_SOURCE,
        "normalization": {
            "P_scale_C_m2": p_scale,
            "P_scale_over_eps0_Vpm": e_scale,
        },
        "cases": cases,
        "same_rms_comparisons": comparisons,
        "summary": {
            "n_cases": len(cases),
            "targets": list(TARGET_RMS_N_ND),
            "spectra_per_target": len(ORIENTATION_SPECTRA),
            "maximum_rms_target_error": max(
                abs(c["orientation_metrics"]["rms_n_ND"] - c["target_rms_n_ND"])
                for c in cases
            ),
            "maximum_field_rms_range_over_mean": max(
                c["normalized_field_rms_spread"]["range_over_mean"] for c in comparisons
            ),
            "maximum_field_p99_range_over_mean": max(
                c["normalized_field_p99_spread"]["range_over_mean"] for c in comparisons
            ),
            "maximum_field_max_range_over_mean": max(
                c["normalized_field_max_spread"]["range_over_mean"] for c in comparisons
            ),
            "all_gauss_residuals_below_1e-7": bool(all(c["gauss_relative_residual"] < 1e-7 for c in cases)),
            "tdgl_coupled": False,
            "switching_enabled": False,
        },
        "interpretation_boundary": [
            "All cases share an RMS normal-projection value constrained by v0.1.18, but their harmonic spectra are project-defined geometry hypotheses.",
            "If field RMS is nearly invariant while p99/max fields vary, RMS n_ND is adequate for average-field scaling but insufficient for hotspot prediction.",
            "If field RMS also varies appreciably at fixed RMS n_ND, higher-order morphology descriptors are required even for bulk local-field magnitude.",
            "Absolute fields remain proportional to the normalized frozen beta/OAF polarization amplitudes and are not phase-resolved measured remanent fields.",
        ],
    }

    outdir = ROOT / "outputs" / "v0119_same_rms_spectrum"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report["summary"], indent=2))
    print(json.dumps(comparisons, indent=2))
    print(f"saved: {path}")


if __name__ == "__main__":
    main()
