import pytest
import numpy as np
from tbm_gfs.green_functions.graphene.single_integral import (
    green_function as single_integral_gf,
)
from tbm_gfs.green_functions.graphene.double_integral import (
    green_function as double_integral_gf,
)
from tbm_gfs.constants import (
    FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY,
)


@pytest.fixture
def sample_energies():
    """
    avoid divergencies at 0 and +t, -t
    """
    t0_abs = np.abs(FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY)
    incr = 0.25
    e_max = 3.0
    return np.concatenate(
        (
            np.arange(-e_max, -t0_abs, incr), 
            np.arange(-t0_abs+incr, 0.0, incr), 
            np.arange(incr, t0_abs-incr, incr), 
            np.arange(t0_abs+incr, e_max, incr), 
        )
    )


def test_g00_same_sublattice(
        sample_energies,
        max_relative_integration_comparison_error,
        ):
    for energy in sample_energies:
        single_integral_result = single_integral_gf(energy, m=0, n=0, s1=1, s2=1)
        double_integral_result = double_integral_gf(energy, m=0, n=0, s1=1, s2=1)
        real_part_err = single_integral_result.real / double_integral_result.real
        imag_part_err = single_integral_result.imag / double_integral_result.imag
        print(
            energy,
            single_integral_result,
            double_integral_result,
            real_part_err,
            imag_part_err,
        )
        assert single_integral_result.real == pytest.approx(
            double_integral_result.real, rel=max_relative_integration_comparison_error
        )
        assert single_integral_result.imag == pytest.approx(
            double_integral_result.imag, rel=max_relative_integration_comparison_error
        )


def test_g00_opposite_sublattice(sample_energies, max_relative_integration_comparison_error):
    for energy in sample_energies:
        single_integral_result = single_integral_gf(energy, m=0, n=0, s1=1, s2=2)
        double_integral_result = double_integral_gf(energy, m=0, n=0, s1=1, s2=2)
        real_part_err = single_integral_result.real / double_integral_result.real
        imag_part_err = single_integral_result.imag / double_integral_result.imag
        print(
            energy,
            single_integral_result,
            double_integral_result,
            real_part_err,
            imag_part_err,
        )
        assert single_integral_result.real == pytest.approx(
            double_integral_result.real, rel=max_relative_integration_comparison_error
        )
        assert single_integral_result.imag == pytest.approx(
            double_integral_result.imag, rel=max_relative_integration_comparison_error
        )



def test_g10_same_sublattice(sample_energies, max_relative_integration_comparison_error):
    for energy in sample_energies:
        single_integral_result = single_integral_gf(energy, m=1, n=0, s1=1, s2=1)
        double_integral_result = double_integral_gf(energy, m=1, n=0, s1=1, s2=1)
        real_part_err = single_integral_result.real / double_integral_result.real
        imag_part_err = single_integral_result.imag / double_integral_result.imag
        print(
            energy,
            single_integral_result,
            double_integral_result,
            real_part_err,
            imag_part_err,
        )
        assert single_integral_result.real == pytest.approx(
            double_integral_result.real, rel=max_relative_integration_comparison_error
        )
        assert single_integral_result.imag == pytest.approx(
            double_integral_result.imag, rel=max_relative_integration_comparison_error
        )


def test_g10_opposite_sublattice(sample_energies, max_relative_integration_comparison_error):
    for energy in sample_energies:
        single_integral_result = single_integral_gf(energy, m=1, n=0, s1=1, s2=2)
        double_integral_result = double_integral_gf(energy, m=1, n=0, s1=1, s2=2)
        real_part_err = single_integral_result.real / double_integral_result.real
        imag_part_err = single_integral_result.imag / double_integral_result.imag
        print(
            energy,
            single_integral_result,
            double_integral_result,
            real_part_err,
            imag_part_err,
        )
        assert single_integral_result.real == pytest.approx(
            double_integral_result.real, rel=max_relative_integration_comparison_error
        )
        assert single_integral_result.imag == pytest.approx(
            double_integral_result.imag, rel=max_relative_integration_comparison_error
        )
