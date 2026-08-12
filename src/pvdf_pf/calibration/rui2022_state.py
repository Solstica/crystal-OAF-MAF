"""Temperature-dependent source states from digitized Rui et al. 2022 SI curves.

This module turns the traceable Figure S1-S3 digitization into temperature-queryable
constitutive targets. It deliberately preserves the source's scope:

* n(T), m_d(T), g(T), lambda(T), and mu_r(T) are film/amorphous-model quantities;
  they are NOT promoted to MOAF-local quantities.
* epsilon_MOAF(T) is a source-model-derived MOAF-specific dielectric target.
* x_RAF(T)/x_MAF(T), required to reconstruct ROAF/MOAF volume fractions, come from
  Figure 5B of the main article and are not contained in the uploaded SI.
* rotational mobility mu_r is not a Debye relaxation time and is not converted to
  tau without an independently validated constitutive relation.

Interpolation is restricted to the digitized temperature range. Linear interpolation
is used for ordinary axes; mu_r is interpolated in log10 space because Figure S3B is
plotted on a logarithmic ordinate and spans several decades.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from pvdf_pf.calibration.source_data import DigitizedPoint, load_digitized_curve_csv


@dataclass(frozen=True)
class Rui2022TemperatureState:
    sample_state: str
    temperature_C: float
    mobility_boundary: str
    active_dipole_density_m3: float
    dipole_moment_D: float
    kirkwood_g: float
    epsilon_moaf: float
    interaction_lambda: float
    rotational_mobility_m2_per_Vs: float
    active_dipole_devitrification_index: float
    rotational_mobility_ratio_from_minus30C: float
    rotational_mobility_decades_from_minus30C: float

    def to_dict(self) -> dict:
        out = asdict(self)
        out["source_scope"] = {
            "epsilon_moaf": "MOAF-specific source-model-derived dielectric target",
            "n_m_d_g": "film/amorphous-model quantities; not MOAF-local",
            "lambda_mu_r": "BOPVDF model quantities; not MOAF-local",
            "devitrification_index": (
                "diagnostic normalized from n(T); not a ROAF/MOAF volume fraction"
            ),
            "mu_r_to_tau": "disabled: no validated conversion is supplied by this module",
        }
        return out


def _filter_points(
    points: list[DigitizedPoint],
    *,
    quantity: str,
    sample_state: str,
    panel: str | None = None,
) -> list[DigitizedPoint]:
    selected = [
        p
        for p in points
        if p.quantity == quantity
        and p.sample_state.lower() == sample_state.lower()
        and (panel is None or p.panel == panel)
    ]
    selected.sort(key=lambda p: p.temperature_C)
    if len(selected) < 2:
        qualifier = f", panel={panel!r}" if panel is not None else ""
        raise ValueError(
            f"need >=2 source points for quantity={quantity!r}, "
            f"sample_state={sample_state!r}{qualifier}"
        )
    temperatures = np.asarray([p.temperature_C for p in selected], dtype=float)
    if np.unique(temperatures).size != temperatures.size:
        raise ValueError("duplicate source temperatures in one curve")
    return selected


def interpolate_digitized_curve(
    points: list[DigitizedPoint],
    temperature_C: float,
    *,
    quantity: str,
    sample_state: str,
    panel: str | None = None,
    log_value: bool = False,
) -> float:
    """Interpolate one digitized source curve without extrapolation."""
    selected = _filter_points(
        points, quantity=quantity, sample_state=sample_state, panel=panel
    )
    t = np.asarray([p.temperature_C for p in selected], dtype=float)
    y = np.asarray([p.value for p in selected], dtype=float)
    target = float(temperature_C)
    if not np.isfinite(target):
        raise ValueError("temperature_C must be finite")
    tol = 1e-12
    if target < t[0] - tol or target > t[-1] + tol:
        raise ValueError(
            f"temperature {target:g} C is outside source range [{t[0]:g}, {t[-1]:g}] C"
        )
    target = min(max(target, float(t[0])), float(t[-1]))
    if log_value:
        if np.any(y <= 0.0):
            raise ValueError("log interpolation requires positive values")
        return float(10.0 ** np.interp(target, t, np.log10(y)))
    return float(np.interp(target, t, y))


def _curve_endpoints(
    points: list[DigitizedPoint], *, quantity: str, sample_state: str, panel: str | None = None
) -> tuple[float, float]:
    selected = _filter_points(
        points, quantity=quantity, sample_state=sample_state, panel=panel
    )
    return float(selected[0].value), float(selected[-1].value)


def build_rui2022_temperature_state(
    data_dir: str | Path,
    *,
    temperature_C: float,
    sample_state: str,
    mobility_boundary: str = "l1",
) -> Rui2022TemperatureState:
    """Build a source-constrained state at one temperature.

    Parameters
    ----------
    data_dir:
        Directory containing figure_s1_digitized.csv, figure_s2_digitized.csv,
        and figure_s3_digitized.csv.
    sample_state:
        ``unpoled`` or ``poled``.
    mobility_boundary:
        ``l1`` (stick) or ``l2`` (slip), matching Figure S3B labels.
    """
    sample = str(sample_state).strip().lower()
    if sample not in {"unpoled", "poled"}:
        raise ValueError("sample_state must be 'unpoled' or 'poled'")
    boundary = str(mobility_boundary).strip().lower()
    if boundary not in {"l1", "l2"}:
        raise ValueError("mobility_boundary must be 'l1' or 'l2'")

    root = Path(data_dir)
    s1 = load_digitized_curve_csv(root / "figure_s1_digitized.csv")
    s2 = load_digitized_curve_csv(root / "figure_s2_digitized.csv")
    s3 = load_digitized_curve_csv(root / "figure_s3_digitized.csv")
    T = float(temperature_C)

    n_1e27 = interpolate_digitized_curve(
        s1, T, quantity="n", sample_state=sample, panel="A"
    )
    m_d = interpolate_digitized_curve(
        s1, T, quantity="m_d", sample_state=sample, panel="B"
    )
    g = interpolate_digitized_curve(
        s1, T, quantity="g", sample_state=sample, panel="C"
    )
    eps_moaf = interpolate_digitized_curve(
        s2, T, quantity="epsilon_MOAF", sample_state=sample, panel="main"
    )
    interaction_lambda = interpolate_digitized_curve(
        s3, T, quantity="lambda", sample_state=sample, panel="A"
    )
    mu_panel = "B_l1" if boundary == "l1" else "B_l2"
    mu_r = interpolate_digitized_curve(
        s3,
        T,
        quantity="mu_r",
        sample_state=sample,
        panel=mu_panel,
        log_value=True,
    )

    n_lo, n_hi = _curve_endpoints(s1, quantity="n", sample_state=sample, panel="A")
    if np.isclose(n_hi, n_lo):
        devitrification_index = 0.0
    else:
        devitrification_index = float((n_1e27 - n_lo) / (n_hi - n_lo))
    devitrification_index = float(np.clip(devitrification_index, 0.0, 1.0))

    mu_lo, _ = _curve_endpoints(
        s3, quantity="mu_r", sample_state=sample, panel=mu_panel
    )
    if mu_lo <= 0.0:
        raise ValueError("source mobility baseline must be positive")
    mobility_ratio = float(mu_r / mu_lo)
    mobility_decades = float(np.log10(mobility_ratio))

    return Rui2022TemperatureState(
        sample_state=sample,
        temperature_C=T,
        mobility_boundary=boundary,
        active_dipole_density_m3=float(n_1e27 * 1.0e27),
        dipole_moment_D=m_d,
        kirkwood_g=g,
        epsilon_moaf=eps_moaf,
        interaction_lambda=interaction_lambda,
        rotational_mobility_m2_per_Vs=mu_r,
        active_dipole_devitrification_index=devitrification_index,
        rotational_mobility_ratio_from_minus30C=mobility_ratio,
        rotational_mobility_decades_from_minus30C=mobility_decades,
    )
