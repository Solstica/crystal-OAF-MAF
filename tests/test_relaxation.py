import numpy as np

from pvdf_pf.calibration.relaxation import (
    debye_delta_epsilon,
    infer_debye_strength_from_real_epsilon,
    laminate_complex_permittivity,
    oaf_complex_permittivity,
)


def test_debye_zero_frequency_limit():
    value = debye_delta_epsilon(0.0, delta_eps=12.0, tau_s=0.01)
    assert np.isclose(value.real, 12.0)
    assert np.isclose(value.imag, 0.0)


def test_debye_high_frequency_real_part_decreases():
    low = oaf_complex_permittivity(1.0, eps_fast=3.0, delta_eps=20.0, tau_s=0.01)
    high = oaf_complex_permittivity(1e6, eps_fast=3.0, delta_eps=20.0, tau_s=0.01)
    assert low.real > high.real
    assert np.isclose(high.real, 3.0, atol=1e-6)


def test_laminate_parallel_and_normal_bounds():
    fractions = {"crystal": 0.5, "oaf": 0.3, "maf": 0.2}
    eps = {"crystal": 3.0 + 0j, "oaf": 30.0 + 0j, "maf": 10.0 + 0j}
    ep, en = laminate_complex_permittivity(fractions, eps)
    assert ep.real > en.real
    assert en.real > 0.0


def test_conditional_debye_inversion_recovers_target():
    fractions = {"crystal": 0.52, "oaf": 0.28, "maf": 0.20}
    strength = infer_debye_strength_from_real_epsilon(
        target_eps_real=19.5,
        frequency_hz=10.0,
        tau_s=0.001,
        phase_fractions=fractions,
        crystal_eps=3.0,
        maf_eps=21.5,
        oaf_eps_fast=3.0,
        orientation="parallel",
    )
    oaf = oaf_complex_permittivity(10.0, eps_fast=3.0, delta_eps=strength, tau_s=0.001)
    ep, _ = laminate_complex_permittivity(
        fractions,
        {"crystal": 3.0 + 0j, "oaf": complex(oaf), "maf": 21.5 + 0j},
    )
    assert np.isclose(ep.real, 19.5, atol=1e-8)
