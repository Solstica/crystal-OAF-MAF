"""Source-model utilities for the ROAF/MOAF split introduced in Rui et al. 2022 SI.

The supporting information writes the BOPVDF dielectric response as a four-component
mixture of crystal, rigid OAF (ROAF), mobile OAF (MOAF), and isotropic amorphous
fraction (IAF).  This module encodes only that accounting identity and its algebraic
inversions.  It does not assign temperature-dependent phase fractions by itself.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class FourFractionState:
    crystal: float
    roaf: float
    moaf: float
    iaf: float

    def as_dict(self) -> dict[str, float]:
        return {
            "crystal": float(self.crystal),
            "roaf": float(self.roaf),
            "moaf": float(self.moaf),
            "iaf": float(self.iaf),
        }

    def validate(self, *, atol: float = 1e-10) -> None:
        vals = np.asarray([self.crystal, self.roaf, self.moaf, self.iaf], dtype=float)
        if np.any(~np.isfinite(vals)) or np.any(vals < 0.0):
            raise ValueError("phase fractions must be finite and non-negative")
        if not np.isclose(vals.sum(), 1.0, atol=atol):
            raise ValueError("crystal + ROAF + MOAF + IAF fractions must sum to one")


def split_total_oaf(
    *,
    crystal_fraction: float,
    total_oaf_fraction: float,
    iaf_fraction: float,
    mobile_fraction_within_oaf: float,
) -> FourFractionState:
    """Split total OAF into rigid and mobile parts without changing total OAF.

    `mobile_fraction_within_oaf` must be independently digitized/fitted from source
    information or another calibrated model.  The function never invents its value.
    """
    fm = float(mobile_fraction_within_oaf)
    if not (0.0 <= fm <= 1.0):
        raise ValueError("mobile_fraction_within_oaf must lie in [0, 1]")
    fc = float(crystal_fraction)
    fo = float(total_oaf_fraction)
    fi = float(iaf_fraction)
    if min(fc, fo, fi) < 0.0:
        raise ValueError("base fractions must be non-negative")
    if not np.isclose(fc + fo + fi, 1.0, atol=1e-10):
        raise ValueError("crystal + total OAF + IAF fractions must sum to one")
    state = FourFractionState(
        crystal=fc,
        roaf=fo * (1.0 - fm),
        moaf=fo * fm,
        iaf=fi,
    )
    state.validate()
    return state


def four_component_parallel_epsilon(
    fractions: FourFractionState,
    *,
    eps_crystal: float,
    eps_roaf: float,
    eps_moaf: float,
    eps_iaf: float,
) -> float:
    """Return the linear mixture used by the SI three-phase/OAF-subsplit model."""
    fractions.validate()
    eps = np.asarray([eps_crystal, eps_roaf, eps_moaf, eps_iaf], dtype=float)
    if np.any(~np.isfinite(eps)) or np.any(eps <= 0.0):
        raise ValueError("all component permittivities must be finite and positive")
    f = np.asarray(
        [fractions.crystal, fractions.roaf, fractions.moaf, fractions.iaf], dtype=float
    )
    return float(np.dot(f, eps))


def infer_moaf_epsilon(
    *,
    film_epsilon: float,
    fractions: FourFractionState,
    eps_crystal: float,
    eps_roaf: float,
    eps_iaf: float,
) -> float:
    """Invert the source mixture equation for epsilon_MOAF.

    This inversion is undefined when the MOAF fraction is zero.  Returned values are
    `SOURCE_MODEL_DERIVED`, not direct local measurements.
    """
    fractions.validate()
    if fractions.moaf <= 0.0:
        raise ValueError("MOAF fraction must be positive for epsilon_MOAF inversion")
    vals = np.asarray([film_epsilon, eps_crystal, eps_roaf, eps_iaf], dtype=float)
    if np.any(~np.isfinite(vals)) or np.any(vals <= 0.0):
        raise ValueError("film and known component permittivities must be positive")
    numerator = (
        float(film_epsilon)
        - fractions.crystal * float(eps_crystal)
        - fractions.roaf * float(eps_roaf)
        - fractions.iaf * float(eps_iaf)
    )
    inferred = numerator / fractions.moaf
    if inferred <= 0.0:
        raise ValueError("source-model inversion produced non-positive epsilon_MOAF")
    return float(inferred)
