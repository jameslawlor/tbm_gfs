import numpy as np
import pytest

from tbm_gfs.recursive_methods.matrices.armchair_nanoribbon import (
    unit_cell_hamiltonian as armchair_nanoribbon_unit_cell_hamiltonian,
    LR_connection_matrix as armchair_nanoribbon_connection_matrix_LR,
    RL_connection_matrix as armchair_nanoribbon_connection_matrix_RL,
)
from tbm_gfs.constants import (
    FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY,
    ARMCHAIR_GNR_PI_BOND_ADJUSTMENT_FACTOR,
)

t = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY
t_edge = t * ARMCHAIR_GNR_PI_BOND_ADJUSTMENT_FACTOR


def test_armchair_nanoribbon_unit_cell_hamiltonian():
    # 4-ZGNR, P.87 Stephen's thesis
    expected_matrix = np.array(
        [
            [0, t, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [t, 0, t, 0, 0, 0, 0, t, 0, 0, 0, 0],
            [0, t, 0, t, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, t, 0, t, 0, 0, 0, 0, t, 0, 0],
            [0, 0, 0, t, 0, t, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, t, 0, 0, 0, 0, 0, 0, t_edge],
            [0, 0, 0, 0, 0, 0, 0, t, 0, 0, 0, 0],
            [0, t, 0, 0, 0, 0, t, 0, t, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, t, 0, t, 0, 0],
            [0, 0, 0, t, 0, 0, 0, 0, t, 0, t, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, t, 0, t],
            [0, 0, 0, 0, 0, t_edge, 0, 0, 0, 0, t, 0],
        ]
    )

    result = armchair_nanoribbon_unit_cell_hamiltonian(size=6)

    assert np.array_equal(expected_matrix, result)


def test_armchair_nanoribbon_connection_matrix_LR():
    # 6-AGNR, P.86 Stephen's thesis
    size = 6
    expected_matrix = np.zeros((2 * size, 2 * size))
    expected_matrix[6][0] = t_edge
    expected_matrix[8][2] = t
    expected_matrix[10][4] = t

    result = armchair_nanoribbon_connection_matrix_LR(size=6)

    assert np.array_equal(expected_matrix, result)


def test_armchair_nanoribbon_connection_matrix_RL():
    # 6-AGNR, P.86 Stephen's thesis
    size = 6

    expected_matrix = np.zeros((2 * size, 2 * size))
    expected_matrix[0][6] = t_edge
    expected_matrix[2][8] = t
    expected_matrix[4][10] = t

    result = armchair_nanoribbon_connection_matrix_RL(size=6)

    assert np.array_equal(expected_matrix, result)


@pytest.mark.parametrize("size", [6, 20])
def test_armchair_nanoribbon_connection_matrices_transpose(size):
    m_RL = armchair_nanoribbon_connection_matrix_RL(size)
    m_LR = armchair_nanoribbon_connection_matrix_LR(size)
    assert np.array_equal(m_RL, np.transpose(m_LR))
