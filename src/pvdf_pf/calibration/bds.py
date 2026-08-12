"""Broadband dielectric spectroscopy utilities for v0.1.2.

The module deliberately separates three levels:

1. spectrum-level Debye parameters (eps_inf, delta_eps, tau);
2. Kirkwood-Frohlich combinations that require an independently supplied active
   dipole number density;
3. the paper-specific n(T), m_d(T), g(T), lambda(T), and rotational mobility,
   which are NOT reconstructed unless numerical source data are supplied.

Complex-permittivity convention used here is eps* = eps' - i eps'', so measured
loss is positive and equals -Im(eps*).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy.optimize import least_squares

EPS0 = 8.8541878128e-12
KB = 1.380649e-23


@dataclass(frozen=True)
class DebyeFitResult:
    eps_inf: float
    delta_eps: float
    tau_s: float
    rmse_real: float
    rmse_loss: float
    n_points: int
    success: bool
    message: str

    def to_dict(self) -> dict:
        return asdict(self)


def debye_complex_permittivity(
    frequency_hz: np.ndarray | float,
    *,
    eps_inf: float,
    delta_eps: float,
    tau_s: float,
) -> np.ndarray:
    """Return eps*(f)=eps_inf+delta_eps/(1+i*2*pi*f*tau)."""
    if eps_inf <= 0.0:
        raise ValueError("eps_inf must be positive")
    if delta_eps < 0.0:
        raise ValueError("delta_eps must be non-negative")
    if tau_s <= 0.0:
        raise ValueError("tau_s must be positive")
    f = np.asarray(frequency_hz, dtype=float)
    if np.any(f < 0.0):
        raise ValueError("frequency must be non-negative")
    return eps_inf + delta_eps / (1.0 + 1j * 2.0 * np.pi * f * tau_s)


def _validate_spectrum(
    frequency_hz: np.ndarray,
    epsilon_real: np.ndarray,
    epsilon_loss: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    f = np.asarray(frequency_hz, dtype=float).reshape(-1)
    er = np.asarray(epsilon_real, dtype=float).reshape(-1)
    el = np.asarray(epsilon_loss, dtype=float).reshape(-1)
    if not (f.size == er.size == el.size):
        raise ValueError("frequency, epsilon_real, and epsilon_loss must have the same length")
    if f.size < 4:
        raise ValueError("at least four frequency points are required")
    if np.any(~np.isfinite(f)) or np.any(~np.isfinite(er)) or np.any(~np.isfinite(el)):
        raise ValueError("spectrum contains non-finite values")
    if np.any(f <= 0.0):
        raise ValueError("BDS fit requires strictly positive frequencies")
    if np.any(er <= 0.0) or np.any(el < 0.0):
        raise ValueError("epsilon_real must be positive and epsilon_loss non-negative")
    if np.unique(f).size < 4:
        raise ValueError("at least four unique frequencies are required")
    order = np.argsort(f)
    return f[order], er[order], el[order]


def fit_debye_spectrum(
    frequency_hz: np.ndarray,
    epsilon_real: np.ndarray,
    epsilon_loss: np.ndarray,
    *,
    eps_inf_bounds: tuple[float, float] = (1.0, 200.0),
    delta_eps_bounds: tuple[float, float] = (1e-6, 1e5),
    tau_bounds_s: tuple[float, float] = (1e-12, 1e3),
) -> DebyeFitResult:
    """Fit one Debye relaxation to epsilon' and positive epsilon'' simultaneously.

    This is a spectrum-level fit, not a claim that RAF/OAF dynamics are exactly
    single-Debye.  It provides a controlled first bridge to the time-domain
    auxiliary-polarization equation.
    """
    f, er, el = _validate_spectrum(frequency_hz, epsilon_real, epsilon_loss)

    e_lo, e_hi = map(float, eps_inf_bounds)
    d_lo, d_hi = map(float, delta_eps_bounds)
    t_lo, t_hi = map(float, tau_bounds_s)
    if not (0.0 < e_lo < e_hi and 0.0 < d_lo < d_hi and 0.0 < t_lo < t_hi):
        raise ValueError("invalid fit bounds")

    eps0_guess = float(np.clip(np.min(er), e_lo * 1.001, e_hi / 1.001))
    delta_guess = float(np.clip(np.max(er) - np.min(er), d_lo * 10.0, d_hi / 10.0))
    if np.max(el) > 0.0:
        f_peak = float(f[np.argmax(el)])
        tau_guess = 1.0 / (2.0 * np.pi * f_peak)
    else:
        tau_guess = 1.0 / (2.0 * np.pi * np.sqrt(f.min() * f.max()))
    tau_guess = float(np.clip(tau_guess, t_lo * 10.0, t_hi / 10.0))

    # Optimize positive parameters in log space except eps_inf.
    x0 = np.array([eps0_guess, np.log(delta_guess), np.log(tau_guess)], dtype=float)
    lower = np.array([e_lo, np.log(d_lo), np.log(t_lo)], dtype=float)
    upper = np.array([e_hi, np.log(d_hi), np.log(t_hi)], dtype=float)

    scale_real = max(float(np.ptp(er)), float(np.mean(er)), 1.0)
    scale_loss = max(float(np.max(el)), 0.05 * scale_real, 1e-3)

    def residual(x: np.ndarray) -> np.ndarray:
        eps_inf = x[0]
        delta_eps = float(np.exp(x[1]))
        tau_s = float(np.exp(x[2]))
        pred = debye_complex_permittivity(
            f, eps_inf=eps_inf, delta_eps=delta_eps, tau_s=tau_s
        )
        return np.concatenate(
            [
                (pred.real - er) / scale_real,
                ((-pred.imag) - el) / scale_loss,
            ]
        )

    fit = least_squares(
        residual,
        x0,
        bounds=(lower, upper),
        method="trf",
        x_scale="jac",
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
        max_nfev=10000,
    )
    eps_inf = float(fit.x[0])
    delta_eps = float(np.exp(fit.x[1]))
    tau_s = float(np.exp(fit.x[2]))
    pred = debye_complex_permittivity(f, eps_inf=eps_inf, delta_eps=delta_eps, tau_s=tau_s)
    return DebyeFitResult(
        eps_inf=eps_inf,
        delta_eps=delta_eps,
        tau_s=tau_s,
        rmse_real=float(np.sqrt(np.mean((pred.real - er) ** 2))),
        rmse_loss=float(np.sqrt(np.mean(((-pred.imag) - el) ** 2))),
        n_points=int(f.size),
        success=bool(fit.success),
        message=str(fit.message),
    )


def kirkwood_frohlich_g_mu2(
    *,
    epsilon_static: float,
    epsilon_fast: float,
    active_dipole_number_density_m3: float,
    temperature_K: float,
) -> float:
    r"""Return the Kirkwood-Frohlich product g*mu^2 in SI units (C^2 m^2).

    Uses

      ((eps_s-eps_inf)(2 eps_s+eps_inf)) /
      (eps_s (eps_inf+2)^2)
        = N g mu^2 / (9 eps0 k_B T).

    `active_dipole_number_density_m3` must come from an independent model or
    measurement.  The function does not infer n(T) from a single spectrum.
    """
    es = float(epsilon_static)
    ei = float(epsilon_fast)
    number_density = float(active_dipole_number_density_m3)
    T = float(temperature_K)
    if not (es > ei > 0.0):
        raise ValueError("require epsilon_static > epsilon_fast > 0")
    if number_density <= 0.0 or T <= 0.0:
        raise ValueError("number density and temperature must be positive")
    lhs = ((es - ei) * (2.0 * es + ei)) / (es * (ei + 2.0) ** 2)
    return float(9.0 * EPS0 * KB * T * lhs / number_density)


def kirkwood_g_from_mu(
    *,
    epsilon_static: float,
    epsilon_fast: float,
    active_dipole_number_density_m3: float,
    temperature_K: float,
    dipole_moment_Cm: float,
) -> float:
    """Infer g when an independent molecular/effective dipole moment is supplied."""
    mu = float(dipole_moment_Cm)
    if mu <= 0.0:
        raise ValueError("dipole_moment_Cm must be positive")
    g_mu2 = kirkwood_frohlich_g_mu2(
        epsilon_static=epsilon_static,
        epsilon_fast=epsilon_fast,
        active_dipole_number_density_m3=active_dipole_number_density_m3,
        temperature_K=temperature_K,
    )
    return float(g_mu2 / mu**2)


def effective_dipole_moment_from_g(
    *,
    epsilon_static: float,
    epsilon_fast: float,
    active_dipole_number_density_m3: float,
    temperature_K: float,
    kirkwood_g: float,
) -> float:
    """Infer effective dipole moment when an independent Kirkwood g is supplied."""
    g = float(kirkwood_g)
    if g <= 0.0:
        raise ValueError("kirkwood_g must be positive")
    g_mu2 = kirkwood_frohlich_g_mu2(
        epsilon_static=epsilon_static,
        epsilon_fast=epsilon_fast,
        active_dipole_number_density_m3=active_dipole_number_density_m3,
        temperature_K=temperature_K,
    )
    return float(np.sqrt(g_mu2 / g))
