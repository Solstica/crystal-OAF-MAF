import numpy as np

from pvdf_pf.core.grid import Grid2D
from pvdf_pf.morphology.orientation_constraint import (
    amplitude_for_target_rms_nz,
    sinusoidal_normal_projection_metrics,
)


def test_zero_target_returns_zero_amplitude_and_zero_projection():
    grid = Grid2D(nz=48, nx=48)
    amplitude = amplitude_for_target_rms_nz(grid, 0.0)
    metrics = sinusoidal_normal_projection_metrics(grid, amplitude)
    assert amplitude == 0.0
    assert metrics["rms_nz"] == 0.0
    assert metrics["max_abs_nz"] == 0.0


def test_inverse_geometry_realizes_representative_targets():
    grid = Grid2D(nz=48, nx=48)
    for target in (0.05, 0.10, 0.20, 0.30):
        amplitude = amplitude_for_target_rms_nz(grid, target)
        metrics = sinusoidal_normal_projection_metrics(grid, amplitude)
        assert np.isclose(metrics["rms_nz"], target, rtol=0.0, atol=1e-10)
        assert amplitude > 0.0


def test_equivalent_amplitude_is_monotonic_with_rms_projection():
    grid = Grid2D(nz=48, nx=48)
    targets = (0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20, 0.25, 0.30)
    amplitudes = [amplitude_for_target_rms_nz(grid, t) for t in targets]
    assert all(b > a for a, b in zip(amplitudes[:-1], amplitudes[1:]))
    assert amplitudes[-1] < 0.10


def test_invalid_target_is_rejected():
    grid = Grid2D(nz=48, nx=48)
    for target in (-0.01, 1.0, 1.1):
        try:
            amplitude_for_target_rms_nz(grid, target)
        except ValueError:
            pass
        else:
            raise AssertionError(f"target {target} should have been rejected")
