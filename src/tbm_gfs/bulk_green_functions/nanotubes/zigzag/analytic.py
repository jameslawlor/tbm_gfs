from cmath import pi, cos, sin, acos, exp

from tbm_gfs.bulk_green_functions.functions import (
    Ne_lambda_function,
    phase_lambda_function,
)
from tbm_gfs.constants import ETA, FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY

t0 = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY


def green_function(n_c, E, m, n, s1, s2):
    E += 1j * ETA
    G = 0

    for j in range(1, n_c + 1):
        kz = j * pi / n_c

        q_A = acos(
            (E**2 - t0**2 - 4 * (t0**2) * (cos(kz) ** 2)) / (4 * t0 * t0 * cos(kz))
        )
        if q_A.imag < 0:
            q_A = -q_A

        if (j == n_c / 2.0) and (m + n == 0):
            if s1 == s2:
                Ne = E
            else:
                Ne = t0

            G += (Ne * exp(1j * kz * (m - n))) / (n_c * (E**2 - t0**2))

        else:
            Ne = Ne_lambda_function(E, kz, q_A, m, n, s1, s2)
            phase_term = phase_lambda_function(kz, q_A, m, n)

            G += (1j / (4 * n_c * (t0 * t0))) * (
                (Ne(kz, q_A) * phase_term(kz, q_A, m, n)) / (cos(kz) * sin(q_A))
            )
    return G
