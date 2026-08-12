import numpy as np
import pytest

from pvdf_pf.calibration.oaf_partition import (
    FourFractionState,
    four_component_parallel_epsilon,
    infer_moaf_epsilon,
    split_total_oaf,
)


def test_split_total_oaf_conserves_phase_fraction():
    state = split_total_oaf(
        crystal_fraction=0.52,
        total_oaf_fraction=0.28,
        iaf_fraction=0.20,
        mobile_fraction_within_oaf=0.25,
    )
    assert np.isclose(state.roaf + state.moaf, 0.28)
    assert np.isclose(sum(state.as_dict().values()), 1.0)


def test_moaf_epsilon_inverse_round_trip():
    state = FourFractionState(crystal=0.52, roaf=0.18, moaf=0.10, iaf=0.20)
    film = four_component_parallel_epsilon(
        state,
        eps_crystal=3.0,
        eps_roaf=3.0,
        eps_moaf=27.0,
        eps_iaf=18.0,
    )
    inferred = infer_moaf_epsilon(
        film_epsilon=film,
        fractions=state,
        eps_crystal=3.0,
        eps_roaf=3.0,
        eps_iaf=18.0,
    )
    assert np.isclose(inferred, 27.0, rtol=1e-13)


def test_moaf_inversion_requires_mobile_fraction():
    state = FourFractionState(crystal=0.52, roaf=0.28, moaf=0.0, iaf=0.20)
    with pytest.raises(ValueError):
        infer_moaf_epsilon(
            film_epsilon=12.0,
            fractions=state,
            eps_crystal=3.0,
            eps_roaf=3.0,
            eps_iaf=18.0,
        )
