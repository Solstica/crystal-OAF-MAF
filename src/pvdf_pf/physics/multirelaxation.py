"""Multiple linear dipolar-relaxation channels for semicrystalline PVDF.

v0.1.3 keeps relaxation channels separate from the ferroelectric TDGL order parameter.
Each channel obeys a Debye auxiliary equation and can be masked to a structural
subregion (e.g. MOAF or IAF).  The physical mapping and channel parameters must come
from source data; this module contains no PVDF parameter guesses.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pvdf_pf.physics.dipolar import EPS0, advance_debye_polarization


@dataclass(frozen=True)
class RelaxationChannel:
    name: str
    tau_s: float
    delta_eps: float

    def validate(self) -> None:
        if not self.name:
            raise ValueError("relaxation channel name must not be empty")
        if not np.isfinite(self.tau_s) or self.tau_s <= 0.0:
            raise ValueError("tau_s must be finite and positive")
        if not np.isfinite(self.delta_eps) or self.delta_eps < 0.0:
            raise ValueError("delta_eps must be finite and non-negative")


def advance_relaxation_bank(
    states: dict[str, np.ndarray],
    E_local: np.ndarray,
    *,
    dt_s: float,
    channels: list[RelaxationChannel],
    masks: dict[str, np.ndarray],
    eps0: float = EPS0,
) -> dict[str, np.ndarray]:
    """Advance all named channels and return a new state dictionary."""
    E = np.asarray(E_local, dtype=float)
    out: dict[str, np.ndarray] = {}
    seen: set[str] = set()
    for channel in channels:
        channel.validate()
        if channel.name in seen:
            raise ValueError(f"duplicate relaxation channel {channel.name!r}")
        seen.add(channel.name)
        if channel.name not in states:
            raise KeyError(f"missing state for channel {channel.name!r}")
        if channel.name not in masks:
            raise KeyError(f"missing mask for channel {channel.name!r}")
        state = np.asarray(states[channel.name], dtype=float)
        mask = np.asarray(masks[channel.name], dtype=bool)
        if state.shape != E.shape or mask.shape != E.shape:
            raise ValueError("channel state, mask, and E_local must have identical shapes")
        out[channel.name] = advance_debye_polarization(
            state,
            E,
            dt_s=dt_s,
            tau_s=channel.tau_s,
            delta_eps=channel.delta_eps,
            eps0=eps0,
            mask=mask,
        )
    return out


def total_relaxational_polarization(states: dict[str, np.ndarray]) -> np.ndarray:
    """Sum channel polarizations without adding background or ferroelectric terms."""
    if not states:
        raise ValueError("states must not be empty")
    arrays = [np.asarray(v, dtype=float) for v in states.values()]
    shape = arrays[0].shape
    if any(a.shape != shape for a in arrays):
        raise ValueError("all channel states must have the same shape")
    total = np.zeros(shape, dtype=float)
    for array in arrays:
        total += array
    return total
