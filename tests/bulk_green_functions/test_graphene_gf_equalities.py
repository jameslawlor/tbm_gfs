import pytest
import numpy as np
from tbm_gfs.bulk_green_functions.graphene.double_integral import (
    green_function as double_integral_gf,
)
from tbm_gfs.bulk_green_functions.graphene.single_integral import (
    green_function as single_integral_gf,
)
from tbm_gfs.bulk_green_functions.graphene.stationary_phase_approximation import (
    spa_same_sublattice_armchair_gf as spa_armchair_gf,
)
from tbm_gfs.constants import (
    FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY,
)


t0_abs = np.abs(FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY)
input_energy_max = 3.0 * t0_abs + t0_abs / 2


@pytest.mark.parametrize(
    "input_energy",
    np.array(
        [
            input_energy_max,  # outside spectral band
            3 * t0_abs / 2,  # E> larger than LDOS peak at plus/minus |t|
            -t0_abs / 2,  # E within spectral band, and below Fermi energy
        ]
    ),
)
@pytest.mark.parametrize("m_and_n_vectors", [(0, 0), (1, 1)])
@pytest.mark.parametrize("s1_and_s2_sublattices", [(1, 1), (1, 2), (2, 1)])
def test_single_and_double_gf_similarity(
    input_energy,
    m_and_n_vectors,
    s1_and_s2_sublattices,
    max_relative_integration_comparison_error,
):
    (m, n) = m_and_n_vectors
    (s1, s2) = s1_and_s2_sublattices
    single_integral_result = single_integral_gf(input_energy, m, n, s1, s2)
    double_integral_result = double_integral_gf(input_energy, m, n, s1, s2)
    assert single_integral_result.real == pytest.approx(
        double_integral_result.real, rel=max_relative_integration_comparison_error
    )
    assert single_integral_result.imag == pytest.approx(
        double_integral_result.imag, rel=max_relative_integration_comparison_error
    )


@pytest.mark.parametrize(
    "input_energy", np.linspace(start=-input_energy_max, stop=input_energy_max, num=10)
)
@pytest.mark.parametrize("m_and_n_vectors", [(0, 0), (1, 1), (1, 0), (8, 3)])
@pytest.mark.parametrize("s1_and_s2_sublattices", [(1, 1), (1, 2), (2, 1)])
def test_single_integral_gf_ka_vs_kz_similarity(
    input_energy,
    m_and_n_vectors,
    s1_and_s2_sublattices,
    max_relative_integration_comparison_error,
):
    (m, n) = m_and_n_vectors
    (s1, s2) = s1_and_s2_sublattices
    kz_integral_result = single_integral_gf(
        input_energy, m, n, s1, s2, integration_variable="kz"
    )
    ka_integral_result = single_integral_gf(
        input_energy, m, n, s1, s2, integration_variable="ka"
    )
    assert kz_integral_result.real == pytest.approx(
        ka_integral_result.real, rel=max_relative_integration_comparison_error
    )
    assert kz_integral_result.imag == pytest.approx(
        ka_integral_result.imag, rel=max_relative_integration_comparison_error
    )


@pytest.mark.parametrize(
    "input_energy",
    np.concatenate(
        (  # Avoid divergencies at +-|t| and +-|3t|
            np.linspace(start=-2.8 * t0_abs, stop=-1.1 * t0_abs, num=10),
            np.linspace(start=-0.8 * t0_abs, stop=0.8 * t0_abs, num=10),
            np.linspace(start=1.1 * t0_abs, stop=2.8 * t0_abs, num=10),
        )
    ),
)
@pytest.mark.parametrize(
    "m_and_n_vectors",
    [
        (5, 5),
        (10, 10),
        (20, 20),
    ],
)  # must be armchair direction i.e. m=n
@pytest.mark.parametrize(
    "s1_and_s2_sublattices",
    [
        (1, 1),
    ],
)
def test_stationary_phase_approximation_armchair(
    input_energy,
    m_and_n_vectors,
    s1_and_s2_sublattices,
):
    max_allowed_relative_err = 0.005
    max_allowed_absolute_err = 0.005

    (m, n) = m_and_n_vectors
    (s1, s2) = s1_and_s2_sublattices
    gf_integral_result = single_integral_gf(
        input_energy, m, n, s1, s2, integration_variable="kz"
    )
    spa_result = spa_armchair_gf(input_energy, m, n, s1, s2)
    assert gf_integral_result.real == pytest.approx(
        spa_result.real,
        rel=max_allowed_relative_err,
        abs=max_allowed_absolute_err,
    )
    assert gf_integral_result.imag == pytest.approx(
        spa_result.imag, rel=max_allowed_relative_err, abs=max_allowed_absolute_err
    )
