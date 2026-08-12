"""Time-domain generalized-Debye auxiliary polarization for v0.1.14.

A broad Cole-Cole response is represented by a finite positive bank of ordinary
Debye modes. Each mode obeys

    tau_j dP_j/dt + P_j = eps0 Delta_eps_j E_local.

The exact exponential zero-order-hold update is used for every mode. The bank is
kept separate from the dimensionless ferroelectric TDGL order parameter until a
physical TDGL time scale is independently calibrated.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pvdf_pf.physics.dipolar import EPS0, advance_debye_polarization


@dataclass(frozen=True)
class GeneralizedDebyeBank:
    epsilon_infinity: float
    tau_modes_s: tuple[float, ...]
    delta_epsilon_modes: tuple[float, ...]
    label: str = "combined_OAF_plus_IAF_amorphous_response"

    def validate(self) -> None:
        tau = np.asarray(self.tau_modes_s, dtype=float)
        strength = np.asarray(self.delta_epsilon_modes, dtype=float)
        if not np.isfinite(self.epsilon_infinity) or self.epsilon_infinity <= 0.0:
            raise ValueError("epsilon_infinity must be finite and positive")
        if tau.size == 0 or tau.size != strength.size:
            raise ValueError("tau_modes_s and delta_epsilon_modes must have equal non-zero length")
        if np.any(~np.isfinite(tau)) or np.any(tau <= 0.0):
            raise ValueError("all mode relaxation times must be finite and positive")
        if np.any(~np.isfinite(strength)) or np.any(strength < 0.0):
            raise ValueError("all mode dielectric strengths must be finite and non-negative")

    @property
    def n_modes(self) -> int:
        return len(self.tau_modes_s)

    @property
    def delta_epsilon_total(self) -> float:
        return float(sum(self.delta_epsilon_modes))


def continuous_generalized_susceptibility(
    bank: GeneralizedDebyeBank,
    frequency_Hz: np.ndarray | float,
) -> np.ndarray:
    bank.validate()
    f = np.asarray(frequency_Hz, dtype=float)
    if np.any(f <= 0.0):
        raise ValueError("frequency must be positive")
    tau = np.asarray(bank.tau_modes_s, dtype=float)
    strength = np.asarray(bank.delta_epsilon_modes, dtype=float)
    omega = 2.0 * np.pi * f
    return np.sum(
        strength / (1.0 + 1j * omega[..., None] * tau),
        axis=-1,
    )


def continuous_generalized_permittivity(
    bank: GeneralizedDebyeBank,
    frequency_Hz: np.ndarray | float,
) -> np.ndarray:
    return float(bank.epsilon_infinity) + continuous_generalized_susceptibility(bank, frequency_Hz)


def advance_generalized_debye_modes(
    P_modes: np.ndarray,
    E_local: np.ndarray | float,
    *,
    dt_s: float,
    bank: GeneralizedDebyeBank,
    mask: np.ndarray | None = None,
    eps0: float = EPS0,
) -> np.ndarray:
    """Advance all auxiliary modes with one exact zero-order-hold step.

    `P_modes` must have shape `(n_modes, *E_local.shape)`. For a scalar field it
    therefore has shape `(n_modes,)`. A spatial mask is broadcast over all modes.
    """
    bank.validate()
    E = np.asarray(E_local, dtype=float)
    P = np.asarray(P_modes, dtype=float)
    expected = (bank.n_modes,) + E.shape
    if P.shape != expected:
        raise ValueError(f"P_modes shape must be {expected}, got {P.shape}")
    tau = np.asarray(bank.tau_modes_s, dtype=float).reshape(
        (bank.n_modes,) + (1,) * E.ndim
    )
    strength = np.asarray(bank.delta_epsilon_modes, dtype=float).reshape(
        (bank.n_modes,) + (1,) * E.ndim
    )
    E_bank = np.broadcast_to(E, expected)
    mode_mask = None
    if mask is not None:
        m = np.asarray(mask, dtype=bool)
        if m.shape != E.shape:
            raise ValueError("mask shape must match E_local")
        mode_mask = np.broadcast_to(m, expected)
    return advance_debye_polarization(
        P,
        E_bank,
        dt_s=dt_s,
        tau_s=tau,
        delta_eps=strength,
        eps0=eps0,
        mask=mode_mask,
    )


def total_auxiliary_polarization(P_modes: np.ndarray) -> np.ndarray:
    P = np.asarray(P_modes, dtype=float)
    if P.ndim < 1 or P.shape[0] == 0:
        raise ValueError("P_modes must contain at least one mode")
    return np.sum(P, axis=0)


def discrete_generalized_susceptibility(
    bank: GeneralizedDebyeBank,
    *,
    frequency_Hz: float,
    dt_s: float,
) -> complex:
    """Exact harmonic susceptibility of the zero-order-hold mode bank."""
    bank.validate()
    f = float(frequency_Hz)
    dt = float(dt_s)
    if f <= 0.0 or dt <= 0.0:
        raise ValueError("frequency_Hz and dt_s must be positive")
    tau = np.asarray(bank.tau_modes_s, dtype=float)
    strength = np.asarray(bank.delta_epsilon_modes, dtype=float)
    a = np.exp(-dt / tau)
    theta = 2.0 * np.pi * f * dt
    return complex(np.sum(strength * (1.0 - a) / (1.0 - a * np.exp(-1j * theta))))
