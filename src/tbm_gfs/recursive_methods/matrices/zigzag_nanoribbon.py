import numpy as np

from tbm_gfs.constants import FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY

t = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY


def unit_cell_hamiltonian(size: int):
    """
    See P.89 in Stephen's thesis

    Note: simpler unit cell structure
    than armchair configuration

    """

    H_uc = np.zeros(shape=(2 * size, 2 * size))

    for ix in range(2 * size - 1):
        H_uc[ix][ix + 1] = t

    H_uc = H_uc + np.transpose(H_uc)

    return H_uc


def LR_connection_matrix(size):
    matrix = np.zeros(shape=(2 * size, 2 * size))

    # black -> white
    for ix in range(0, 2 * size, 4):
        matrix[ix][ix + 1] = t

    # white -> black
    for ix in range(0, 2 * size, 4):
        matrix[ix + 3][ix + 2] = t

    return matrix


def RL_connection_matrix(size):
    matrix = np.zeros(shape=(2 * size, 2 * size))

    # white -> black
    for ix in range(0, 2 * size, 4):
        matrix[ix + 1][ix] = t

    # black -> white
    for ix in range(0, 2 * size, 4):
        matrix[ix + 2][ix + 3] = t

    return matrix
