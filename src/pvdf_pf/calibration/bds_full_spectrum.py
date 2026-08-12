"""v0.1.13 same-state BOPVDF full-spectrum relaxation fitting.

Figure 2 alpha-window digitized spectra are fitted at the FILM level first.
The measured/static film permittivity is fixed from Figure 3A. Two models are
compared: Debye and Cole-Cole. A small non-negative loss floor is treated as a
nuisance term only; it is not promoted to a material parameter.

Only after a film-level fit is accepted is the parallel two-phase relation used
to map the fit to the combined amorphous response (OAF+IAF). No OAF-only
response is inferred here.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import csv

import numpy as np
from scipy.optimize import least_squares

ETA_CRYSTAL = 0.52
EPSILON_CRYSTAL = 3.0


@dataclass(frozen=True)
class Figure2SpectrumPoint:
    sample_state: str
    temperature_C: float
    frequency_Hz: float
    epsilon_real: float
    epsilon_loss: float
    digitization_score_real: float
    digitization_score_loss: float
    provenance: str


@dataclass(frozen=True)
class RelaxationFit:
    model: str
    sample_state: str
    temperature_C: float
    n_points: int
    epsilon_static_film: float
    epsilon_infinity_film: float
    delta_epsilon_film: float
    tau_s: float
    beta_cole_cole: float
    alpha_cole_cole: float
    loss_floor: float
    rmse_real: float
    rmse_loss: float
    aic: float
    raw_peak_frequency_Hz: float
    peak_coverage_ratio: float
    coverage_status: str
    epsilon_static_amorphous: float
    epsilon_infinity_amorphous: float
    delta_epsilon_amorphous: float
    response_assignment: str = "combined_OAF_plus_IAF_amorphous_response"

    def to_dict(self) -> dict:
        return asdict(self)


def load_figure2_alpha_window(path: str | Path) -> list[Figure2SpectrumPoint]:
    p = Path(path)
    out: list[Figure2SpectrumPoint] = []
    with p.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        required = {
            "sample_state", "temperature_C", "frequency_Hz", "epsilon_real",
            "epsilon_loss", "digitization_score_real", "digitization_score_loss",
            "selection", "provenance",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing Figure-2 spectrum columns: {sorted(missing)}")
        for row in reader:
            if row["selection"] != "ALPHA_WINDOW":
                continue
            vals = [float(row[k]) for k in (
                "temperature_C", "frequency_Hz", "epsilon_real", "epsilon_loss",
                "digitization_score_real", "digitization_score_loss"
            )]
            if not np.all(np.isfinite(vals)):
                raise ValueError("non-finite Figure-2 digitized value")
            if vals[1] <= 0 or vals[2] <= 0 or vals[3] < 0:
                raise ValueError("invalid Figure-2 spectrum value")
            if row["provenance"] != "DIGITIZED_SOURCE":
                raise ValueError("Figure-2 plot transcription must retain DIGITIZED_SOURCE provenance")
            out.append(Figure2SpectrumPoint(
                sample_state=row["sample_state"].strip().lower(),
                temperature_C=vals[0], frequency_Hz=vals[1],
                epsilon_real=vals[2], epsilon_loss=vals[3],
                digitization_score_real=vals[4], digitization_score_loss=vals[5],
                provenance=row["provenance"],
            ))
    if not out:
        raise ValueError("no Figure-2 alpha-window points found")
    return out


def cole_cole_complex(
    frequency_Hz,
    *,
    epsilon_static: float,
    epsilon_infinity: float,
    tau_s: float,
    beta: float,
) -> np.ndarray:
    """Return eps*=eps_inf+(eps_s-eps_inf)/(1+(i*w*tau)^beta).

    beta=1 is Debye. The common Cole-Cole alpha is alpha_CC=1-beta.
    Convention: eps*=eps'-i eps''.
    """
    f = np.asarray(frequency_Hz, dtype=float)
    es = float(epsilon_static)
    ei = float(epsilon_infinity)
    tau = float(tau_s)
    b = float(beta)
    if np.any(f <= 0):
        raise ValueError("frequency must be positive")
    if not (es > ei > 0):
        raise ValueError("require epsilon_static > epsilon_infinity > 0")
    if tau <= 0 or not (0 < b <= 1):
        raise ValueError("invalid tau/beta")
    w = 2.0 * np.pi * f
    return ei + (es - ei) / (1.0 + (1j * w * tau) ** b)


def _combined_amorphous(
    value_film: float,
    *,
    eta_crystal: float = ETA_CRYSTAL,
    epsilon_crystal: float = EPSILON_CRYSTAL,
) -> float:
    return float((value_film - eta_crystal * epsilon_crystal) / (1.0 - eta_crystal))


def fit_relaxation(
    points: list[Figure2SpectrumPoint],
    *,
    epsilon_static_film: float,
    model: str = "cole-cole",
    max_digitization_score: float = 90.0,
) -> RelaxationFit:
    if model not in {"debye", "cole-cole"}:
        raise ValueError("model must be 'debye' or 'cole-cole'")
    kept = [
        p for p in points
        if p.digitization_score_real < max_digitization_score
        and p.digitization_score_loss < max_digitization_score
    ]
    if len(kept) < 8:
        raise ValueError("insufficient quality-controlled Figure-2 points")
    kept.sort(key=lambda p: p.frequency_Hz)
    sample = kept[0].sample_state
    temperature = kept[0].temperature_C
    if any(p.sample_state != sample or p.temperature_C != temperature for p in kept):
        raise ValueError("fit_relaxation requires one sample_state/temperature group")

    f = np.asarray([p.frequency_Hz for p in kept])
    er = np.asarray([p.epsilon_real for p in kept])
    el = np.asarray([p.epsilon_loss for p in kept])
    es = float(epsilon_static_film)
    if es <= 1:
        raise ValueError("epsilon_static_film must be >1")

    raw_peak = float(f[np.argmax(el)])
    epsinf0 = float(np.clip(np.median(er[-min(5, len(er)):]), 1.01, es - 0.05))
    tau0 = 1.0 / (2.0 * np.pi * raw_peak)
    beta_free = model == "cole-cole"
    if beta_free:
        x0 = np.array([epsinf0, np.log(tau0), 0.45, max(0.0, float(np.min(el)) * 0.2)])
        lo = np.array([1.0, np.log(1e-12), 0.15, 0.0])
        hi = np.array([es - 1e-4, np.log(1e2), 1.0, max(5.0, float(np.max(el)))])
    else:
        x0 = np.array([epsinf0, np.log(tau0), max(0.0, float(np.min(el)) * 0.2)])
        lo = np.array([1.0, np.log(1e-12), 0.0])
        hi = np.array([es - 1e-4, np.log(1e2), max(5.0, float(np.max(el)))])

    scale_r = max(float(np.ptp(er)), 1.0)
    scale_l = max(float(np.max(el)), 0.3)

    def residual(x):
        ei = float(x[0])
        tau = float(np.exp(x[1]))
        beta = float(x[2]) if beta_free else 1.0
        floor = float(x[3] if beta_free else x[2])
        pred = cole_cole_complex(
            f,
            epsilon_static=es,
            epsilon_infinity=ei,
            tau_s=tau,
            beta=beta,
        )
        return np.r_[
            (pred.real - er) / scale_r,
            ((-pred.imag + floor) - el) / scale_l,
        ]

    opt = least_squares(
        residual,
        x0,
        bounds=(lo, hi),
        method="trf",
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
        max_nfev=20000,
    )
    ei = float(opt.x[0])
    tau = float(np.exp(opt.x[1]))
    beta = float(opt.x[2]) if beta_free else 1.0
    floor = float(opt.x[3] if beta_free else opt.x[2])
    pred = cole_cole_complex(
        f,
        epsilon_static=es,
        epsilon_infinity=ei,
        tau_s=tau,
        beta=beta,
    )
    pred_loss = -pred.imag + floor
    rmse_r = float(np.sqrt(np.mean((pred.real - er) ** 2)))
    rmse_l = float(np.sqrt(np.mean((pred_loss - el) ** 2)))
    sse = float(np.sum((pred.real - er) ** 2) + np.sum((pred_loss - el) ** 2))
    nobs = 2 * len(f)
    k = len(opt.x)
    aic = float(nobs * np.log(max(sse / nobs, 1e-30)) + 2 * k)

    coverage = float(1.0e7 / raw_peak)
    status = "RIGHT_CENSORED_PEAK" if coverage < 3.0 else "PEAK_COVERAGE_ACCEPTABLE"
    es_am = _combined_amorphous(es)
    ei_am = _combined_amorphous(ei)

    return RelaxationFit(
        model=model,
        sample_state=sample,
        temperature_C=temperature,
        n_points=len(f),
        epsilon_static_film=es,
        epsilon_infinity_film=ei,
        delta_epsilon_film=float(es - ei),
        tau_s=tau,
        beta_cole_cole=beta,
        alpha_cole_cole=float(1.0 - beta),
        loss_floor=floor,
        rmse_real=rmse_r,
        rmse_loss=rmse_l,
        aic=aic,
        raw_peak_frequency_Hz=raw_peak,
        peak_coverage_ratio=coverage,
        coverage_status=status,
        epsilon_static_amorphous=es_am,
        epsilon_infinity_amorphous=ei_am,
        delta_epsilon_amorphous=float(es_am - ei_am),
    )
