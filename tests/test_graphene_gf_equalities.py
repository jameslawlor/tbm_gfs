import pytest
import numpy as np
from tbm_gfs.green_functions.graphene.single_integral_kz import (
    green_function as single_integral_gf_kz,
)
from tbm_gfs.green_functions.graphene.single_integral_ka import (
    green_function as single_integral_gf_ka,
)
from tbm_gfs.green_functions.graphene.double_integral import (
    green_function as double_integral_gf,
)
from tbm_gfs.constants import (
    FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY,
)


t0_abs = np.abs(FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY)
e_max = 3.0 + t0_abs / 2


@pytest.mark.parametrize("input_energy", np.linspace(start=-e_max, stop=e_max, num=10))
@pytest.mark.parametrize("m_and_n_vectors", [(0, 0), (1, 1), (1, 0)])
@pytest.mark.parametrize("s1_and_s2_sublattices", [(1, 1), (1, 2), (2, 1)])
def test_single_and_double_gf_similarity(
    input_energy,
    m_and_n_vectors,
    s1_and_s2_sublattices,
    max_relative_integration_comparison_error,
):
    (m, n) = m_and_n_vectors
    (s1, s2) = s1_and_s2_sublattices
    single_integral_result = single_integral_gf_kz(input_energy, m, n, s1, s2)
    double_integral_result = double_integral_gf(input_energy, m, n, s1, s2)
    assert single_integral_result.real == pytest.approx(
        double_integral_result.real, rel=max_relative_integration_comparison_error
    )
    assert single_integral_result.imag == pytest.approx(
        double_integral_result.imag, rel=max_relative_integration_comparison_error
    )


@pytest.mark.parametrize("input_energy", np.linspace(start=-e_max, stop=e_max, num=10))
@pytest.mark.parametrize("m_and_n_vectors", [(0, 0), (1, 1), (1, 0)])
@pytest.mark.parametrize("s1_and_s2_sublattices", [(1, 1), (1, 2), (2, 1)])
def test_single_integral_gf_ka_vs_kz_similarity(
    input_energy,
    m_and_n_vectors,
    s1_and_s2_sublattices,
    max_relative_integration_comparison_error,
):
    (m, n) = m_and_n_vectors
    (s1, s2) = s1_and_s2_sublattices
    kz_integral_result = single_integral_gf_kz(input_energy, m, n, s1, s2)
    ka_integral_result = single_integral_gf_ka(input_energy, m, n, s1, s2)
    assert kz_integral_result.real == pytest.approx(
        ka_integral_result.real, rel=max_relative_integration_comparison_error
    )
    assert kz_integral_result.imag == pytest.approx(
        ka_integral_result.imag, rel=max_relative_integration_comparison_error
    )
