import numpy as np
import pytest

from tbm_gfs.recursive_methods.matrices.zigzag_nanoribbon import (
    unit_cell_hamiltonian as zigzag_nanoribbon_unit_cell_hamiltonian,
    LR_connection_matrix as zigzag_nanoribbon_connection_matrix_LR,
    RL_connection_matrix as zigzag_nanoribbon_connection_matrix_RL,
)
from tbm_gfs.constants import (
    FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY,
)

t = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY


def test_zigzag_nanoribbon_unit_cell_hamiltonian():
    # 4-ZGNR, P.87 Stephen's thesis
    expected_matrix = np.array(
        [
            [
                0,
                t,
                0,
                0,
                0,
                0,
                0,
                0,
            ],
            [
                t,
                0,
                t,
                0,
                0,
                0,
                0,
                0,
            ],
            [
                0,
                t,
                0,
                t,
                0,
                0,
                0,
                0,
            ],
            [
                0,
                0,
                t,
                0,
                t,
                0,
                0,
                0,
            ],
            [
                0,
                0,
                0,
                t,
                0,
                t,
                0,
                0,
            ],
            [
                0,
                0,
                0,
                0,
                t,
                0,
                t,
                0,
            ],
            [
                0,
                0,
                0,
                0,
                0,
                t,
                0,
                t,
            ],
            [
                0,
                0,
                0,
                0,
                0,
                0,
                t,
                0,
            ],
        ]
    )

    result = zigzag_nanoribbon_unit_cell_hamiltonian(size=4)

    assert np.array_equal(expected_matrix, result)


def test_zigzag_nanoribbon_connection_matrix_LR():
    size = 4
    expected_matrix = np.zeros((2 * size, 2 * size))
    expected_matrix[0][1] = t
    expected_matrix[3][2] = t
    expected_matrix[4][5] = t
    expected_matrix[7][6] = t

    result = zigzag_nanoribbon_connection_matrix_LR(size=size)

    assert np.array_equal(expected_matrix, result)


def test_zigzag_nanoribbon_connection_matrix_RL():
    size = 4

    expected_matrix = np.zeros((2 * size, 2 * size))
    expected_matrix[1][0] = t
    expected_matrix[2][3] = t
    expected_matrix[5][4] = t
    expected_matrix[6][7] = t
    result = zigzag_nanoribbon_connection_matrix_RL(size=size)

    assert np.array_equal(expected_matrix, result)


@pytest.mark.parametrize("size", [6, 20])
def test_zigzag_nanoribbon_connection_matrices_transpose(size):
    m_RL = zigzag_nanoribbon_connection_matrix_RL(size)
    m_LR = zigzag_nanoribbon_connection_matrix_LR(size)
    assert np.array_equal(m_RL, np.transpose(m_LR))
