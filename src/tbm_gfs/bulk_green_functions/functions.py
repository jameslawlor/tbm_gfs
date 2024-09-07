from typing import Callable
from cmath import cos, exp, sqrt
import numpy as np

from tbm_gfs.constants import (
    FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY,
)

t0 = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY

# TODO: Refactor lambda functions to clean up unneeded arguments


def Ne_lambda_function(
    E: complex, kz: complex, ka: complex, m: int, n: int, s1: int, s2: int
) -> Callable[[complex, complex], complex]:
    if s1 == s2:
        return lambda kz, ka: E
    elif s1 == 1 and s2 == 2:
        return lambda kz, ka: t0 * (1 + 2 * cos(kz) * exp(1j * ka))
    elif s1 == 2 and s2 == 1:
        exponent_sign = 1 if (m == 0 and n == 0) else -1
        return lambda kz, ka: t0 * (1 + 2 * cos(kz) * exp(exponent_sign * 1j * ka))
    else:
        raise ValueError("Error! Check inputs s1 and s2 are valid. Must equal 1 or 2. ")


def phase_lambda_function(
    kz: complex,
    ka: complex,
    m: int,
    n: int,
) -> Callable[[complex, complex, int, int], complex]:
    return lambda kz, ka, m, n: exp(1j * (ka * (m + n) + kz * (m - n)))


def dispersion_relation_lambda_function(
    kz: complex, ka: complex
) -> Callable[[complex, complex], complex]:
    return lambda kz, ka: t0 * sqrt(1 + 4 * cos(ka) * cos(kz) + 4 * cos(kz) * cos(kz))


def graphene_sector_method(m, n, s1, s2):
    """
    Transform any vector in the graphene lattice
    into the irreducible sector from 0 to pi/6

    Parameters:
    - m, n: Integers representing the coordinates in the lattice.
    - s1, s2: Integers representing the sublattice indices (0 or 1).

    Returns:
    - A tuple (m, n, s1, s2) where (m, n) lies in the irreducible sector.
    """

    s = s2 - s1

    if s < 0:  # white-black case
        m, n, s = -n, -m, -s

    if m >= 0:
        if n >= 0:  # sector 0
            pass

        elif abs(n) <= m:  # sector 1
            m, n, s = -n, m + n + s, -s

        elif abs(n) > m:  # sector 2
            m, n, s = m, -n - m + -s, s
    elif m < 0:
        if n <= 0:  # sector 3
            m, n, s = -n, -m, -s

        elif abs(m) > n:  # sector 4
            m, n, s = n, -n - m - s, s

        elif abs(m) <= n:  # sector 5
            m, n, s = -m, n + m + (s), -s
    if n > m:  # to irreducible sector
        m, n = n, m

    if s == 0:
        s1 = 1
        s2 = 1
    elif s > 0:
        s1 = 1
        s2 = 2
    elif s < 0:
        s1 = 2
        s2 = 1

    return m, n, s1, s2


def convert_lattice_vectors_to_real_space(m, n, s):
    """
    Units of nearest-neighbour distance a=1
    """

    a1 = [3 / 2, np.sqrt(3) / 2]
    a2 = [3 / 2, -np.sqrt(3) / 2]

    x_coordinate = m * a1[0] + n * a2[0]
    y_coordinate = m * a1[1] + n * a2[1]

    if s == 2:
        x_coordinate += 1

    return x_coordinate, y_coordinate
