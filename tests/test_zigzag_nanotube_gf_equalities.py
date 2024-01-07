import pytest
import numpy as np
from tbm_gfs.green_functions.nanotubes.zigzag.analytic import (
    green_function as zigzag_analytic_gf,
)
from tbm_gfs.green_functions.nanotubes.zigzag.integral import (
    green_function as zigzag_integral_gf,
)
from tbm_gfs.constants import (
    FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY,
)


t0_abs = np.abs(FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY)


@pytest.mark.parametrize(
    "input_energy", np.linspace(start=-3.0 * t0_abs, stop=3.0 * t0_abs, num=13)
)
@pytest.mark.parametrize(
    "nanotube_cirumferences",
    [
        6,
        11,
    ],
)
@pytest.mark.parametrize(
    "m_and_n_vectors",
    [
        (0, 0),
        (1, 1),
        (1, 0),
        (8, 3),
    ],
)
@pytest.mark.parametrize(
    "s1_and_s2_sublattices",
    [
        (1, 1),
        (1, 2),
        (2, 1),
    ],
)
def test_zigzag_nanotube_analytic_vs_integral_gf_similarity(
    input_energy,
    nanotube_cirumferences,
    m_and_n_vectors,
    s1_and_s2_sublattices,
    max_relative_integration_comparison_error,
):
    (m, n) = m_and_n_vectors
    (s1, s2) = s1_and_s2_sublattices
    n_c = nanotube_cirumferences
    zigzag_analytic_result = zigzag_analytic_gf(
        n_c,
        input_energy,
        m,
        n,
        s1,
        s2,
    )
    zigzag_integral_result = zigzag_integral_gf(
        n_c,
        input_energy,
        m,
        n,
        s1,
        s2,
    )
    assert zigzag_analytic_result.real == pytest.approx(
        zigzag_integral_result.real,
        rel=max_relative_integration_comparison_error,
        abs=max_relative_integration_comparison_error,
    )
    assert zigzag_analytic_result.imag == pytest.approx(
        zigzag_integral_result.imag,
        rel=max_relative_integration_comparison_error,
        abs=max_relative_integration_comparison_error,
    )
