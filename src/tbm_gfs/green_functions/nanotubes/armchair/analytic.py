from cmath import pi, cos, sin, acos, sqrt

from tbm_gfs.green_functions.functions import (
    Ne_lambda_function,
    phase_lambda_function,
)
from tbm_gfs.constants import ETA, FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY

t0 = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY


def green_function(n_c, E, m, n, s1, s2):
    E += 1j * ETA

    G = 0

    for j in range(1, n_c + 1):
        ka = pi * j / n_c

        for sign in (-1, 1):
            q_Z = acos(
                (-1 / 2) * cos(ka)
                + (sign / (2 * t0)) * sqrt(E**2 - (t0**2) * sin(ka) ** 2)
            )

            if q_Z.imag < 0:
                q_Z = -q_Z

            Ne = Ne_lambda_function(E, q_Z, ka, m, n, s1, s2)
            phase_term = phase_lambda_function(q_Z, ka, m, n)

            G += (
                Ne(q_Z, ka)
                * (((1j) / (4 * n_c * t0)) * phase_term(q_Z, ka, m, n))
                / (sin(q_Z) * (t0 * cos(ka) + 2 * t0 * cos(q_Z)))
            )
    return G
