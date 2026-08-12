from __future__ import annotations

import numpy as np


def debye_delta_epsilon(
    frequency_hz: np.ndarray | float,
    *,
    delta_eps: float,
    tau_s: float,
) -> np.ndarray:
    r"""Complex Debye contribution Δε/(1 + iωτ).

    This v0.1 utility is a constitutive scaffold.  `delta_eps` and `tau_s` are not
    treated as physical OAF parameters unless a provenance record explicitly marks
    them as fitted or directly supported.
    """
    if delta_eps < 0.0:
        raise ValueError("delta_eps must be non-negative")
    if tau_s <= 0.0:
        raise ValueError("tau_s must be positive")
    f = np.asarray(frequency_hz, dtype=float)
    if np.any(f < 0.0):
        raise ValueError("frequency must be non-negative")
    omega = 2.0 * np.pi * f
    return delta_eps / (1.0 + 1j * omega * tau_s)


def oaf_complex_permittivity(
    frequency_hz: np.ndarray | float,
    *,
    eps_fast: float,
    delta_eps: float,
    tau_s: float,
) -> np.ndarray:
    """Return ε*_OAF(ω) = ε_fast + Debye dipolar contribution."""
    if eps_fast <= 0.0:
        raise ValueError("eps_fast must be positive")
    return eps_fast + debye_delta_epsilon(
        frequency_hz,
        delta_eps=delta_eps,
        tau_s=tau_s,
    )


def laminate_complex_permittivity(
    phase_fractions: dict[str, float],
    phase_eps_complex: dict[str, complex],
) -> tuple[complex, complex]:
    r"""Return ideal-laminate parallel and normal complex permittivities.

    For a perfectly flat laminate,

        ε_parallel = Σ f_i ε_i,
        ε_normal   = 1 / Σ f_i / ε_i.

    These are geometry limits used for identifiability checks, not a substitute for
    the 2-D heterogeneous cell problem.
    """
    names = ("crystal", "oaf", "maf")
    fractions = np.asarray([phase_fractions[name] for name in names], dtype=float)
    if np.any(fractions <= 0.0) or not np.isclose(fractions.sum(), 1.0, atol=1e-10):
        raise ValueError("phase fractions must be positive and sum to one")
    eps = np.asarray([phase_eps_complex[name] for name in names], dtype=complex)
    if np.any(np.real(eps) <= 0.0):
        raise ValueError("real part of each phase permittivity must be positive")

    eps_parallel = complex(np.sum(fractions * eps))
    eps_normal = complex(1.0 / np.sum(fractions / eps))
    return eps_parallel, eps_normal


def infer_debye_strength_from_real_epsilon(
    *,
    target_eps_real: float,
    frequency_hz: float,
    tau_s: float,
    phase_fractions: dict[str, float],
    crystal_eps: float,
    maf_eps: float,
    oaf_eps_fast: float,
    orientation: str = "parallel",
) -> float:
    """Infer Δε_OAF from one real-permittivity target under a chosen τ hypothesis.

    A single frequency cannot identify both Δε and τ.  This function therefore
    requires τ to be supplied and returns the corresponding Δε.  Results from this
    function must be labelled `INFERRED_CONDITIONAL`.
    """
    if target_eps_real <= 0.0:
        raise ValueError("target_eps_real must be positive")
    if tau_s <= 0.0:
        raise ValueError("tau_s must be positive")
    if frequency_hz < 0.0:
        raise ValueError("frequency_hz must be non-negative")
    if orientation not in {"parallel", "normal"}:
        raise ValueError("orientation must be 'parallel' or 'normal'")

    def real_eff(delta_eps: float) -> float:
        eps_oaf = complex(
            oaf_complex_permittivity(
                frequency_hz,
                eps_fast=oaf_eps_fast,
                delta_eps=delta_eps,
                tau_s=tau_s,
            )
        )
        ep, en = laminate_complex_permittivity(
            phase_fractions,
            {"crystal": complex(crystal_eps), "oaf": eps_oaf, "maf": complex(maf_eps)},
        )
        return float(np.real(ep if orientation == "parallel" else en))

    base = real_eff(0.0)
    if target_eps_real < base - 1e-12:
        raise ValueError(
            f"target real permittivity {target_eps_real:g} is below the zero-relaxation baseline {base:g}"
        )

    lo = 0.0
    hi = 1.0
    while real_eff(hi) < target_eps_real and hi < 1.0e6:
        hi *= 2.0
    if hi >= 1.0e6 and real_eff(hi) < target_eps_real:
        raise ValueError("target could not be reached with delta_eps <= 1e6")

    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if real_eff(mid) < target_eps_real:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
