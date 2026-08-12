"""Positive generalized-Debye (Prony) approximation for v0.1.14.

The accepted v0.1.13 Cole-Cole response is a frequency-domain constitutive fit.
This module converts that broad, causal response into a finite bank of ordinary
Debye modes with non-negative dielectric strengths.  The modes are numerical
quadrature/Prony modes, not independently identified molecular mechanisms.

All conversions operate on the *combined OAF+IAF amorphous response* inherited
from v0.1.13.  Nothing in this module promotes the result to an OAF-only model.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy.optimize import nnls

from pvdf_pf.calibration.bds_full_spectrum import RelaxationFit, cole_cole_complex


@dataclass(frozen=True)
class PronyRepresentation:
    sample_state: str
    temperature_C: float
    epsilon_infinity_amorphous: float
    epsilon_static_amorphous: float
    target_tau_s: float
    target_beta_cole_cole: float
    frequency_min_Hz: float
    frequency_max_Hz: float
    n_modes: int
    tau_modes_s: tuple[float, ...]
    delta_epsilon_modes: tuple[float, ...]
    normalized_rms_complex_error: float
    normalized_max_complex_error: float
    static_strength_error: float
    active_mode_count: int
    response_assignment: str = "combined_OAF_plus_IAF_amorphous_response"
    mode_interpretation: str = "positive numerical generalized-Debye modes; not physical subphases"

    @property
    def delta_epsilon_total(self) -> float:
        return float(sum(self.delta_epsilon_modes))

    def to_dict(self) -> dict:
        return asdict(self)


def generalized_debye_susceptibility(
    frequency_Hz: np.ndarray | float,
    *,
    tau_modes_s: np.ndarray | tuple[float, ...],
    delta_epsilon_modes: np.ndarray | tuple[float, ...],
) -> np.ndarray:
    """Return sum_j Delta_eps_j/(1+i*omega*tau_j)."""
    f = np.asarray(frequency_Hz, dtype=float)
    tau = np.asarray(tau_modes_s, dtype=float).reshape(-1)
    strength = np.asarray(delta_epsilon_modes, dtype=float).reshape(-1)
    if tau.size == 0 or tau.size != strength.size:
        raise ValueError("tau_modes_s and delta_epsilon_modes must be non-empty and equal length")
    if np.any(f <= 0.0) or np.any(tau <= 0.0) or np.any(strength < 0.0):
        raise ValueError("frequency/tau must be positive and strengths non-negative")
    omega = 2.0 * np.pi * f
    return np.sum(
        strength / (1.0 + 1j * omega[..., None] * tau),
        axis=-1,
    )


def generalized_debye_permittivity(
    frequency_Hz: np.ndarray | float,
    *,
    epsilon_infinity: float,
    tau_modes_s: np.ndarray | tuple[float, ...],
    delta_epsilon_modes: np.ndarray | tuple[float, ...],
) -> np.ndarray:
    if epsilon_infinity <= 0.0:
        raise ValueError("epsilon_infinity must be positive")
    return float(epsilon_infinity) + generalized_debye_susceptibility(
        frequency_Hz,
        tau_modes_s=tau_modes_s,
        delta_epsilon_modes=delta_epsilon_modes,
    )


def fit_positive_prony_from_cole_cole(
    fit: RelaxationFit,
    *,
    frequency_min_Hz: float = 1e-2,
    frequency_max_Hz: float = 1e9,
    n_modes: int = 49,
    support_margin_decades: float = 1.5,
    n_fit_frequencies: int = 601,
    static_constraint_weight: float = 100.0,
) -> PronyRepresentation:
    """Approximate a v0.1.13 Cole-Cole fit by passive Debye modes.

    The fitting interval brackets the 1--1e7 Hz source measurement window by two
    additional frequency decades in total.  Mode characteristic frequencies are
    extended further by `support_margin_decades` so the finite bank can represent
    the Cole-Cole tails without assigning negative weights.

    Non-negative least squares is followed by an exact normalization of the mode
    weights, so sum(Delta_eps_j) equals the v0.1.13 amorphous dielectric strength.
    """
    if fit.model != "cole-cole":
        raise ValueError("v0.1.14 Prony conversion requires an accepted Cole-Cole fit")
    fmin = float(frequency_min_Hz)
    fmax = float(frequency_max_Hz)
    if not (0.0 < fmin < fmax):
        raise ValueError("require 0 < frequency_min_Hz < frequency_max_Hz")
    if n_modes < 5 or n_fit_frequencies < max(50, n_modes):
        raise ValueError("insufficient Prony modes or fitting frequencies")
    if support_margin_decades < 0.0 or static_constraint_weight <= 0.0:
        raise ValueError("support margin must be non-negative and constraint weight positive")

    delta_total = float(fit.delta_epsilon_amorphous)
    if delta_total <= 0.0:
        raise ValueError("amorphous dielectric strength must be positive")

    f = np.logspace(np.log10(fmin), np.log10(fmax), int(n_fit_frequencies))
    f_modes = np.logspace(
        np.log10(fmin) - float(support_margin_decades),
        np.log10(fmax) + float(support_margin_decades),
        int(n_modes),
    )
    tau_modes = np.sort(1.0 / (2.0 * np.pi * f_modes))

    target = cole_cole_complex(
        f,
        epsilon_static=fit.epsilon_static_amorphous,
        epsilon_infinity=fit.epsilon_infinity_amorphous,
        tau_s=fit.tau_s,
        beta=fit.beta_cole_cole,
    )
    target_norm = (
        target - fit.epsilon_infinity_amorphous
    ) / delta_total

    omega = 2.0 * np.pi * f
    basis = 1.0 / (1.0 + 1j * omega[:, None] * tau_modes[None, :])
    A = np.vstack([
        basis.real,
        -basis.imag,
        float(static_constraint_weight) * np.ones((1, n_modes)),
    ])
    y = np.concatenate([
        target_norm.real,
        -target_norm.imag,
        np.array([float(static_constraint_weight)]),
    ])
    weights, _ = nnls(A, y)
    weight_sum = float(weights.sum())
    if not np.isfinite(weight_sum) or weight_sum <= 0.0:
        raise RuntimeError("positive Prony fit produced zero/non-finite total weight")
    weights /= weight_sum
    strengths = delta_total * weights

    pred_norm = basis @ weights
    err = pred_norm - target_norm
    rms = float(np.sqrt(np.mean(np.abs(err) ** 2)))
    maxerr = float(np.max(np.abs(err)))
    static_error = float(strengths.sum() - delta_total)
    active = int(np.count_nonzero(strengths > delta_total * 1e-8))

    return PronyRepresentation(
        sample_state=fit.sample_state,
        temperature_C=float(fit.temperature_C),
        epsilon_infinity_amorphous=float(fit.epsilon_infinity_amorphous),
        epsilon_static_amorphous=float(fit.epsilon_static_amorphous),
        target_tau_s=float(fit.tau_s),
        target_beta_cole_cole=float(fit.beta_cole_cole),
        frequency_min_Hz=fmin,
        frequency_max_Hz=fmax,
        n_modes=int(n_modes),
        tau_modes_s=tuple(float(v) for v in tau_modes),
        delta_epsilon_modes=tuple(float(v) for v in strengths),
        normalized_rms_complex_error=rms,
        normalized_max_complex_error=maxerr,
        static_strength_error=static_error,
        active_mode_count=active,
    )
