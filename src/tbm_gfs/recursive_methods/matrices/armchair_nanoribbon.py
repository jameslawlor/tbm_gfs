import numpy as np

from tbm_gfs.constants import FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY, ARMCHAIR_GNR_PI_BOND_ADJUSTMENT_FACTOR

t = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY
t_edge = t * ARMCHAIR_GNR_PI_BOND_ADJUSTMENT_FACTOR

def unit_cell_hamiltonian(size: int):
    """
    See P.86 in Stephen's thesis

    Method used here is based on the observation
    that the matrix is composed of two groups of two identical
    smaller sub-matrices.

    """
    # fill top left quadrant, corresponds to connecting the atoms
    # lying across the width of the ribbon
    H_upper_left = np.zeros(shape=(size, size))
    for ix in range(size - 1):
        H_upper_left[ix][ix + 1] = t
    for ix in range(1, size):
        H_upper_left[ix][ix - 1] = t

    H_bottom_right = H_upper_left

    # Now generate top right quadrant
    # corresponds to length-wise atom
    # connections within unit cell
    H_upper_right = np.zeros(shape=(size, size))
    for ix in range(1, size, 2):
        H_upper_right[ix][ix] = t

    H_upper_right[-1][-1] = t_edge

    H_bottom_left = H_upper_right

    H_uc = np.block([[H_upper_left, H_upper_right], [H_bottom_left, H_bottom_right]])

    return H_uc


def LR_connection_matrix(size):
    matrix = np.zeros(shape=(2 * size, 2 * size))
    for ix in range(0, size, 2):
        matrix[size + ix][ix] = t

    matrix[size][0] = t_edge
    return matrix


def RL_connection_matrix(size):
    matrix = np.zeros(shape=(2 * size, 2 * size))

    for ix in range(0, size, 2):
        matrix[ix][size + ix] = t

    matrix[0][size] = t_edge
    return matrix
