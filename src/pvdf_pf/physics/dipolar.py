"""Time-domain auxiliary dipolar polarization for v0.1.2.

This module implements the exact-in-time update of the linear Debye equation

    tau dP_rel/dt + P_rel = eps0 * Delta_eps * E_local.

The auxiliary polarization is kept separate from the existing dimensionless TDGL
order parameter because the TDGL time unit has not yet been calibrated to seconds.
It can be coupled to TDGL only after a physical time-scale mapping is established.
"""

from __future__ import annotations

import numpy as np

EPS0 = 8.8541878128e-12


def advance_debye_polarization(
    P_rel: np.ndarray | float,
    E_local: np.ndarray | float,
    *,
    dt_s: float,
    tau_s: np.ndarray | float,
    delta_eps: np.ndarray | float,
    eps0: float = EPS0,
    mask: np.ndarray | None = None,
) -> np.ndarray:
    """Advance one Debye auxiliary-polarization step using the exact exponential map.

    E_local is assumed constant during the step.  `mask` can restrict the response
    to OAF cells.  Outside the mask the returned auxiliary polarization is zero.
    """
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive")
    if eps0 <= 0.0:
        raise ValueError("eps0 must be positive")

    P = np.asarray(P_rel, dtype=float)
    E = np.asarray(E_local, dtype=float)
    tau = np.asarray(tau_s, dtype=float)
    strength = np.asarray(delta_eps, dtype=float)
    if np.any(tau <= 0.0):
        raise ValueError("tau_s must be positive")
    if np.any(strength < 0.0):
        raise ValueError("delta_eps must be non-negative")

    decay = np.exp(-dt_s / tau)
    target = eps0 * strength * E
    updated = decay * P + (1.0 - decay) * target

    if mask is not None:
        m = np.asarray(mask, dtype=bool)
        try:
            updated = np.where(m, updated, 0.0)
        except ValueError as exc:
            raise ValueError("mask is not broadcast-compatible with polarization") from exc
    return np.asarray(updated, dtype=float)


def simulate_debye_history(
    E_history: np.ndarray,
    *,
    dt_s: float,
    tau_s: float,
    delta_eps: float,
    P0_Cpm2: float = 0.0,
    eps0: float = EPS0,
) -> np.ndarray:
    """Integrate a scalar electric-field history for validation and BDS bridging."""
    field = np.asarray(E_history, dtype=float).reshape(-1)
    if field.size == 0:
        raise ValueError("E_history must not be empty")
    out = np.empty_like(field, dtype=float)
    P = float(P0_Cpm2)
    for i, E in enumerate(field):
        P = float(
            advance_debye_polarization(
                P,
                E,
                dt_s=dt_s,
                tau_s=tau_s,
                delta_eps=delta_eps,
                eps0=eps0,
            )
        )
        out[i] = P
    return out


def complex_susceptibility_from_harmonic_history(
    E_history: np.ndarray,
    P_history: np.ndarray,
    *,
    dt_s: float,
    frequency_hz: float,
    eps0: float = EPS0,
    discard_fraction: float = 0.5,
) -> complex:
    """Estimate Delta-epsilon* from a steady harmonic time-domain trajectory.

    The convention is epsilon* = epsilon' - i epsilon''.  The returned quantity is
    the auxiliary dipolar contribution P/(eps0 E), not total permittivity.
    """
    E = np.asarray(E_history, dtype=float).reshape(-1)
    P = np.asarray(P_history, dtype=float).reshape(-1)
    if E.size != P.size or E.size < 8:
        raise ValueError("E_history and P_history must have equal length >= 8")
    if dt_s <= 0.0 or frequency_hz <= 0.0 or eps0 <= 0.0:
        raise ValueError("dt_s, frequency_hz and eps0 must be positive")
    if not (0.0 <= discard_fraction < 1.0):
        raise ValueError("discard_fraction must lie in [0, 1)")

    start = int(np.floor(discard_fraction * E.size))
    E = E[start:]
    P = P[start:]
    t = dt_s * np.arange(start, start + E.size)
    basis = np.exp(-1j * 2.0 * np.pi * frequency_hz * t)
    Ehat = np.sum(E * basis)
    Phat = np.sum(P * basis)
    if abs(Ehat) < 1e-30:
        raise ValueError("harmonic electric-field amplitude is numerically zero")
    # With e^{-i wt} projection and E=cos(wt), the dynamic equation produces a
    # phasor proportional to 1/(1+i wtau), matching the repository convention.
    return complex(Phat / (eps0 * Ehat))
