# import numpy as np
# import pytest

from tbm_gfs.recursive_methods.matrices.armchair_nanoribbon import (
    unit_cell_hamiltonian as armchair_nanoribbon_unit_cell_hamiltonian,
    LR_connection_matrix as armchair_nanoribbon_connection_matrix_LR,
    RL_connection_matrix as armchair_nanoribbon_connection_matrix_RL,
)

from tbm_gfs.recursive_methods.matrices.zigzag_nanoribbon import (
    unit_cell_hamiltonian as zigzag_nanoribbon_unit_cell_hamiltonian,
    LR_connection_matrix as zigzag_nanoribbon_connection_matrix_LR,
    RL_connection_matrix as zigzag_nanoribbon_connection_matrix_RL,
)
from tbm_gfs.constants import (
    FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY,
)

t = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY