"""Four-component small-signal dielectric cell model for v0.1.7.

This module combines:
- source-model crystal/ROAF permittivity assumptions (3.0),
- digitized/source-model-derived epsilon_MOAF(T),
- epsilon_IAF(T) obtained by the SI-requested melt extrapolation,
- a spatial crystal/ROAF/MOAF/IAF map.

The result is a linear dielectric calibration problem only. These permittivities
contain orientational response and are not TDGL background permittivities.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from pvdf_pf.calibration.dielectric import effective_permittivity
from pvdf_pf.calibration.iaf_permittivity import epsilon_iaf_from_melt_extrapolation
from pvdf_pf.calibration.rui2022_state import build_rui2022_temperature_state
from pvdf_pf.core.spectral import grad_periodic
from pvdf_pf.morphology.oaf_devitrification import (
    SOURCE_CRYSTAL,
    SOURCE_ROAF,
    SOURCE_MOAF,
    SOURCE_IAF,
)


@dataclass(frozen=True)
class FourPhasePermittivityState:
    temperature_C: float
    epsilon_crystal: float
    epsilon_roaf: float
    epsilon_moaf: float
    epsilon_iaf: float
    sample_state: str
    iaf_provenance: str = "SOURCE_MODEL_EXTRAPOLATION"
    moaf_provenance: str = "SOURCE_MODEL_DERIVED"
    crystal_roaf_provenance: str = "SOURCE_MODEL_ASSUMPTION"

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class LocalFieldSummary:
    epsilon_effective: float
    field_mean: float
    field_std: float
    field_p05: float
    field_p95: float
    field_max: float
    crystal_mean: float
    roaf_mean: float
    moaf_mean: float
    iaf_mean: float

    def to_dict(self) -> dict:
        return asdict(self)


def build_four_phase_permittivity_state(
    rui2022_dir: str | Path,
    *,
    temperature_C: float,
    sample_state: str = "poled",
    epsilon_crystal: float = 3.0,
    epsilon_roaf: float = 3.0,
) -> FourPhasePermittivityState:
    root = Path(rui2022_dir)
    source_state = build_rui2022_temperature_state(
        root,
        temperature_C=float(temperature_C),
        sample_state=sample_state,
        mobility_boundary="l1",
    )
    eps_iaf, _ = epsilon_iaf_from_melt_extrapolation(
        root / "figure_1b_melt_eps_digitized.csv",
        temperature_C=float(temperature_C),
    )
    return FourPhasePermittivityState(
        temperature_C=float(temperature_C),
        epsilon_crystal=float(epsilon_crystal),
        epsilon_roaf=float(epsilon_roaf),
        epsilon_moaf=float(source_state.epsilon_moaf),
        epsilon_iaf=float(eps_iaf),
        sample_state=str(sample_state).strip().lower(),
    )


def four_phase_eps_map(
    source_map: np.ndarray,
    state: FourPhasePermittivityState,
) -> np.ndarray:
    phase = np.asarray(source_map)
    allowed = {SOURCE_CRYSTAL, SOURCE_ROAF, SOURCE_MOAF, SOURCE_IAF}
    if not set(np.unique(phase)).issubset(allowed):
        raise ValueError("source_map contains unsupported phase codes")
    out = np.empty(phase.shape, dtype=float)
    out[phase == SOURCE_CRYSTAL] = state.epsilon_crystal
    out[phase == SOURCE_ROAF] = state.epsilon_roaf
    out[phase == SOURCE_MOAF] = state.epsilon_moaf
    out[phase == SOURCE_IAF] = state.epsilon_iaf
    return out


def cell_center_field(
    psi: np.ndarray,
    grid,
    *,
    axis: str,
    E0: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Return cell-centered (Ez, Ex) from the dielectric corrector potential."""
    dpsi_dz, dpsi_dx = grad_periodic(psi, grid.dz, grid.dx)
    Ez = -dpsi_dz
    Ex = -dpsi_dx
    if axis == "z":
        Ez = Ez + float(E0)
    elif axis == "x":
        Ex = Ex + float(E0)
    else:
        raise ValueError("axis must be 'z' or 'x'")
    return Ez, Ex


def _phase_mean(field: np.ndarray, phase: np.ndarray, code: int) -> float:
    mask = phase == code
    return float(np.mean(field[mask])) if np.any(mask) else float("nan")


def solve_four_phase_cell(
    source_map: np.ndarray,
    state: FourPhasePermittivityState,
    grid,
    *,
    axis: str = "z",
    E0: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, LocalFieldSummary]:
    """Solve the linear heterogeneous dielectric cell problem and summarize |E|/E0."""
    eps_map = four_phase_eps_map(source_map, state)
    eps_eff, psi = effective_permittivity(eps_map, grid, axis=axis, E0=E0)
    Ez, Ex = cell_center_field(psi, grid, axis=axis, E0=E0)
    emag = np.sqrt(Ez**2 + Ex**2) / abs(float(E0))
    phase = np.asarray(source_map)
    summary = LocalFieldSummary(
        epsilon_effective=float(eps_eff),
        field_mean=float(np.mean(emag)),
        field_std=float(np.std(emag)),
        field_p05=float(np.quantile(emag, 0.05)),
        field_p95=float(np.quantile(emag, 0.95)),
        field_max=float(np.max(emag)),
        crystal_mean=_phase_mean(emag, phase, SOURCE_CRYSTAL),
        roaf_mean=_phase_mean(emag, phase, SOURCE_ROAF),
        moaf_mean=_phase_mean(emag, phase, SOURCE_MOAF),
        iaf_mean=_phase_mean(emag, phase, SOURCE_IAF),
    )
    return eps_map, Ez, Ex, summary
