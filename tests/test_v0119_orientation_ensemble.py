import numpy as np

from pvdf_pf.core.grid import Grid2D
from pvdf_pf.morphology.orientation_ensemble import (
    HarmonicMode,
    scale_for_target_rms_normal,
    spectral_orientation_metrics,
    spectral_three_phase,
)
from pvdf_pf.morphology.three_phase import CRYSTAL, OAF, MAF


def families():
    return (
        (HarmonicMode(1, 1.0, 0.0),),
        (HarmonicMode(2, 1.0, 0.0),),
        (HarmonicMode(1, 1.0, 0.0), HarmonicMode(2, 0.5, 1.0)),
        (
            HarmonicMode(1, 1.0, 0.2),
            HarmonicMode(2, 0.5, 1.1),
            HarmonicMode(3, 0.3, 2.0),
            HarmonicMode(4, 0.15, 2.7),
        ),
    )


def test_all_spectra_can_realize_same_rms_orientation():
    grid = Grid2D(nz=48, nx=48)
    target = 0.15
    for modes in families():
        scale = scale_for_target_rms_normal(grid, modes, target)
        metrics = spectral_orientation_metrics(grid, modes, scale)
        assert np.isclose(metrics["rms_n_ND"], target, rtol=0.0, atol=1e-10)
        assert abs(metrics["mean_n_ND"]) < 2e-2


def test_spectral_phase_map_contains_all_structural_regions():
    grid = Grid2D(nz=48, nx=48)
    modes = families()[-1]
    scale = scale_for_target_rms_normal(grid, modes, 0.20)
    phase, n_nd = spectral_three_phase(
        grid,
        crystal_fraction=5.78 / 11.8,
        oaf_fraction=3.02 / 11.8,
        iaf_fraction=3.00 / 11.8,
        modes=modes,
        scale=scale,
    )
    assert phase.shape == grid.shape
    assert n_nd.shape == grid.shape
    assert np.any(phase == CRYSTAL)
    assert np.any(phase == OAF)
    assert np.any(phase == MAF)
    assert np.isclose(np.sqrt(np.mean(n_nd**2)), 0.20, atol=1e-10)


def test_higher_mode_requires_smaller_coordinate_scale_for_same_rms():
    grid = Grid2D(nz=48, nx=48)
    target = 0.10
    scale_k1 = scale_for_target_rms_normal(grid, (HarmonicMode(1, 1.0),), target)
    scale_k2 = scale_for_target_rms_normal(grid, (HarmonicMode(2, 1.0),), target)
    assert scale_k2 < scale_k1
    assert np.isclose(scale_k2 / scale_k1, 0.5, rtol=0.03)


def test_invalid_harmonic_mode_is_rejected():
    grid = Grid2D(nz=48, nx=48)
    try:
        scale_for_target_rms_normal(grid, (HarmonicMode(0, 1.0),), 0.10)
    except ValueError:
        pass
    else:
        raise AssertionError("zero harmonic mode should have been rejected")
