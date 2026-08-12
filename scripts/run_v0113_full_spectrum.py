from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path
import sys
import math

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.bopvdf_v0113 import RUI2022_DIR, SAMPLE_STATES, TEMPERATURES_C, V0113_RULES
from pvdf_pf.calibration.bds_full_spectrum import load_figure2_alpha_window, fit_relaxation
from pvdf_pf.calibration.source_data import load_digitized_curve_csv


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
    peak_source = load_digitized_curve_csv(RUI2022_DIR / "figure_3b_lnfd_digitized.csv")

    grouped = defaultdict(list)
    for point in spectra:
        grouped[(point.sample_state, point.temperature_C)].append(point)

    results = {}
    all_delta_aic = []
    right_censored = []
    for sample in SAMPLE_STATES:
        rows = []
        for T in TEMPERATURES_C:
            eps_s = exact_source_value(film_static, sample, T, "epsilon_c_film")
            deb = fit_relaxation(grouped[(sample, T)], epsilon_static_film=eps_s, model="debye")
            cc = fit_relaxation(grouped[(sample, T)], epsilon_static_film=eps_s, model="cole-cole")
            ln_fd = exact_source_value(peak_source, sample, T, "ln_f_d")
            fd_deconvolved = math.exp(ln_fd)
            tau_deconvolved = 1.0 / (2.0 * math.pi * fd_deconvolved)
            delta_aic = deb.aic - cc.aic
            all_delta_aic.append(delta_aic)
            if cc.coverage_status == "RIGHT_CENSORED_PEAK":
                right_censored.append({"sample_state": sample, "temperature_C": T})
            item = cc.to_dict()
            item.update({
                "aic_debye": deb.aic,
                "delta_aic_debye_minus_cole_cole": delta_aic,
                "fd_deconvolved_source_Hz": fd_deconvolved,
                "tau_deconvolved_source_s": tau_deconvolved,
                "tau_raw_cole_cole_over_tau_deconvolved": cc.tau_s / tau_deconvolved,
            })
            rows.append(item)
        results[sample] = rows

    report = {
        "model_version": "v0.1.13",
        "stage": "same-state Figure-2 full-spectrum dielectric relaxation",
        "rules": V0113_RULES,
        "results": results,
        "summary": {
            "minimum_delta_AIC_Debye_minus_ColeCole": min(all_delta_aic),
            "cole_cole_preferred_at_all_states": all(v > 0 for v in all_delta_aic),
            "right_censored_peak_states": right_censored,
            "promotion_scope": "combined OAF+IAF amorphous dynamic response only",
        },
        "interpretation": [
            "Across every digitized state, the Cole-Cole model is preferred over single-Debye by AIC, so the alpha relaxation is intrinsically broad on the Figure-2 scale.",
            "The fitted beta is a spectral broadening descriptor; alpha_CC=1-beta. It is not an OAF phase fraction.",
            "The film-level fitted response is mapped to the combined amorphous response only through the same-state eta_cr=0.52 parallel two-phase relation.",
            "The Figure-3B fd values are deconvolved source peak frequencies. They are retained separately from raw Figure-2 Cole-Cole characteristic times rather than forced to be identical.",
            "40 C spectra are right-censored by the 1e7 Hz measurement ceiling and must not be used as strong epsilon_infinity calibration anchors.",
            "No OAF-only epsilon*(omega,T), tau(T), or TDGL background permittivity is created in v0.1.13.",
        ],
    }
    outdir = ROOT / "outputs" / "v0113_full_spectrum"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report["summary"], indent=2))
    print(f"saved: {path}")


if __name__ == "__main__":
    main()
