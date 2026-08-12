from __future__ import annotations

from collections import defaultdict
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v0114 import PRONY, RUI2022_DIR, SAMPLE_STATES, TEMPERATURES_C, V0114_RULES
from pvdf_pf.calibration.bds_full_spectrum import load_figure2_alpha_window, fit_relaxation
from pvdf_pf.calibration.prony import fit_positive_prony_from_cole_cole
from pvdf_pf.calibration.source_data import load_digitized_curve_csv
from pvdf_pf.physics.generalized_debye import (
    GeneralizedDebyeBank,
    continuous_generalized_susceptibility,
    discrete_generalized_susceptibility,
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


def main() -> None:
    spectra = load_figure2_alpha_window(RUI2022_DIR / "figure_2_bds_alpha_window_digitized.csv")
    film_static = load_digitized_curve_csv(RUI2022_DIR / "figure_3a_film_eps_digitized.csv")
    grouped = defaultdict(list)
    for point in spectra:
        grouped[(point.sample_state, point.temperature_C)].append(point)

    prony_kwargs = {k: v for k, v in PRONY.items() if k != "max_normalized_complex_error"}
    results = {}
    max_rms = 0.0
    max_abs = 0.0
    max_discrete_error = 0.0
    min_active = 10**9
    all_passive = True
    all_static_closed = True

    for sample in SAMPLE_STATES:
        rows = []
        for T in TEMPERATURES_C:
            eps_s = exact_source_value(film_static, sample, T, "epsilon_c_film")
            cc = fit_relaxation(
                grouped[(sample, T)],
                epsilon_static_film=eps_s,
                model="cole-cole",
            )
            rep = fit_positive_prony_from_cole_cole(cc, **prony_kwargs)
            strengths = rep.delta_epsilon_modes
            all_passive = all_passive and all(v >= 0.0 for v in strengths)
            all_static_closed = all_static_closed and abs(rep.static_strength_error) < 1e-10
            max_rms = max(max_rms, rep.normalized_rms_complex_error)
            max_abs = max(max_abs, rep.normalized_max_complex_error)
            min_active = min(min_active, rep.active_mode_count)

            bank = GeneralizedDebyeBank(
                epsilon_infinity=rep.epsilon_infinity_amorphous,
                tau_modes_s=rep.tau_modes_s,
                delta_epsilon_modes=rep.delta_epsilon_modes,
            )
            f_probe = 1.0 / (2.0 * math.pi * rep.target_tau_s)
            dt_s = 1.0 / (2000.0 * f_probe)
            chi_cont = complex(continuous_generalized_susceptibility(bank, f_probe))
            chi_disc = discrete_generalized_susceptibility(
                bank, frequency_Hz=f_probe, dt_s=dt_s
            )
            discrete_error = abs(chi_disc - chi_cont) / rep.delta_epsilon_total
            max_discrete_error = max(max_discrete_error, discrete_error)

            item = rep.to_dict()
            item.update({
                "probe_frequency_Hz": f_probe,
                "probe_dt_s": dt_s,
                "normalized_ZOH_vs_continuous_error_at_probe": discrete_error,
                "coverage_status_inherited_from_v0113": cc.coverage_status,
            })
            rows.append(item)
        results[sample] = rows

    tolerance = float(PRONY["max_normalized_complex_error"])
    report = {
        "model_version": "v0.1.14",
        "stage": "Cole-Cole to positive generalized-Debye time-domain bridge",
        "rules": V0114_RULES,
        "prony_settings": PRONY,
        "results": results,
        "summary": {
            "all_mode_strengths_non_negative": all_passive,
            "static_dielectric_strength_closed": all_static_closed,
            "maximum_normalized_rms_complex_error": max_rms,
            "maximum_normalized_complex_error": max_abs,
            "approximation_tolerance": tolerance,
            "all_states_within_tolerance": max_abs <= tolerance,
            "minimum_active_mode_count": min_active,
            "maximum_normalized_ZOH_vs_continuous_error_at_probe": max_discrete_error,
            "promotion_scope": "combined OAF+IAF amorphous auxiliary polarization only",
        },
        "interpretation": [
            "The broad Cole-Cole response is represented by a passive positive sum of ordinary Debye modes, so no fractional time derivative is required in the local-field solver.",
            "The modal dielectric strengths sum exactly to the v0.1.13 combined-amorphous dielectric strength; the modes are numerical representation modes rather than identified structural subphases.",
            "The exact exponential zero-order-hold update has an explicitly auditable discrete harmonic response; its finite-dt phase/amplitude error is kept separate from the material approximation error.",
            "The representation remains OAF+IAF combined. It is not promoted to an OAF-only response and is not yet coupled to the dimensionless TDGL clock.",
        ],
    }
    outdir = ROOT / "outputs" / "v0114_prony_bridge"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report["summary"], indent=2))
    print(f"saved: {path}")


if __name__ == "__main__":
    main()
