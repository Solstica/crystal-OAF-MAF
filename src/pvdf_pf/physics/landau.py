from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class PhasePolynomial:
    """Scalar local constitutive polynomial used by the v0 numerical baseline."""

    a2: float
    a4: float
    a6: float = 0.0
    eps_b: float = 1.0


def build_phase_parameter_maps(phase_map: np.ndarray, params: dict[str, PhasePolynomial]) -> dict[str, np.ndarray]:
    """Map phase-local coefficients onto the grid. Phase ids: 0 crystal, 1 OAF, 2 MAF."""
    required = ("crystal", "oaf", "maf")
    missing = [name for name in required if name not in params]
    if missing:
        raise KeyError(f"missing phase parameters: {missing}")

    maps = {key: np.zeros(phase_map.shape, dtype=float) for key in ("a2", "a4", "a6", "eps_b")}
    for pid, name in enumerate(required):
        mask = phase_map == pid
        p = params[name]
        maps["a2"][mask] = p.a2
        maps["a4"][mask] = p.a4
        maps["a6"][mask] = p.a6
        maps["eps_b"][mask] = p.eps_b

    if np.any(maps["eps_b"] <= 0.0):
        raise ValueError("background permittivity must remain positive")
    return maps


def local_energy_derivative(P: np.ndarray, maps: dict[str, np.ndarray], E_total_z: np.ndarray | float) -> np.ndarray:
    """d f_local / dP including -E.P; gradient energy is handled separately."""
    return maps["a2"] * P + maps["a4"] * P**3 + maps["a6"] * P**5 - E_total_z


def local_energy_density(P: np.ndarray, maps: dict[str, np.ndarray], E_total_z: np.ndarray | float = 0.0) -> np.ndarray:
    return (
        0.5 * maps["a2"] * P**2
        + 0.25 * maps["a4"] * P**4
        + (maps["a6"] / 6.0) * P**6
        - E_total_z * P
    )
