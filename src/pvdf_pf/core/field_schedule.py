import numpy as np


def triangular_field(n_steps: int, amplitude: float = 1.0, cycles: int = 1) -> np.ndarray:
    """Full -E -> +E -> -E triangular cycles, including both endpoints."""
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")
    if cycles < 1:
        raise ValueError("cycles must be >= 1")
    t = np.linspace(0.0, float(cycles), n_steps + 1, endpoint=True)
    frac = t % 1.0
    tri = np.where(frac < 0.5, -1.0 + 4.0 * frac, 3.0 - 4.0 * frac)
    tri[-1] = -1.0
    return amplitude * tri
