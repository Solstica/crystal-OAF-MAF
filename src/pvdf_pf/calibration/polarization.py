"""Polarization-partition constraints derived from macroscopic BOPVDF data."""

from __future__ import annotations


def required_phase_polarization(
    *,
    film_polarization: float,
    target_fraction: float,
    other_contributions: list[tuple[float, float]],
) -> float:
    """Return the target-phase polarization required by a volume-fraction closure.

    Parameters use SI polarization units consistently, typically C/m^2.
    `other_contributions` contains `(fraction, local_polarization)` pairs.
    """
    if target_fraction <= 0.0:
        raise ValueError("target_fraction must be positive")
    accounted = sum(float(f) * float(p) for f, p in other_contributions)
    return (float(film_polarization) - accounted) / float(target_fraction)


def oaf_polarization_lower_bound(
    *,
    film_ps_Cpm2: float,
    crystal_fraction: float,
    oaf_fraction: float,
    crystal_ps_upper_bound_Cpm2: float,
    maf_ps_Cpm2: float = 0.0,
    maf_fraction: float = 0.0,
) -> float:
    """Lower bound on OAF local Ps under a crystal-Ps upper bound.

    The bound assumes the crystal contributes at its allowed maximum.  Any lower
    crystal polarization would require a larger OAF contribution.
    """
    return required_phase_polarization(
        film_polarization=film_ps_Cpm2,
        target_fraction=oaf_fraction,
        other_contributions=[
            (crystal_fraction, crystal_ps_upper_bound_Cpm2),
            (maf_fraction, maf_ps_Cpm2),
        ],
    )
