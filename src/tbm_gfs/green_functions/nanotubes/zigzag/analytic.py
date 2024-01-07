from cmath import pi, cos, sin, acos, exp

from tbm_gfs.green_functions.functions import (
    Ne_lambda_function,
    phase_lambda_function,
)
from tbm_gfs.constants import ETA, FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY

t0 = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY

def green_function(n_c, E, m, n, S1, S2):
    E += 1j*ETA

    if (m + n) < 0:
        S1, S2 = S2, S1
        m, n = -m, -n
    G = 0

    for j in range(1, n_c + 1):
        kZ = j * pi / n_c

        q = acos((E**2 - t0**2 - 4 * (t0**2) * (cos(kZ) ** 2)) / (4 * t0 * t0 * cos(kZ)))
        if q.imag < 0:
            q = -q

        if (j == n_c / 2.0) and (m + n == 0):
            if S1 == S2:
                Ne = E
            else:
                Ne = t0
            G += (Ne * exp(1j * kZ * (m - n))) / (n_c * (E**2 - t0**2))

        else:
            if S1 == S2:
                Ne = E
            if S1 < S2:
                Ne = t0 + 2 * t0 * cos(kZ) * exp(1j * q)
            if S1 > S2:
                if m == 0 and n == 0:
                    Ne = t0 + 2 * t0 * cos(kZ) * exp(1j * q)
                else:
                    Ne = t0 + 2 * t0 * cos(kZ) * exp(-1j * q)

            G += (1j / (4 * n_c * (t0 * t0))) * (
                (Ne * exp(1j * kZ * (m - n)) * exp(1j * q * (m + n)))
                / (cos(kZ) * sin(q))
            )
    return G
