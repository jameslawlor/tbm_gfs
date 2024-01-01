"""
See Equation 13 in
Duffy, J. M., et al.
"Variable range of the RKKY interaction in edged graphene."
Journal of Physics: Condensed Matter 26.5 (2013): 055007.
"""

from scipy.integrate import quad
from cmath import sin, cos, acos, pi, exp, sqrt
from typing import Callable

from tbm_gfs.constants import (
    FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY,
    ETA,
    INTEGRATION_LIMIT,
    INTEGRATION_EPS_ABS,
    INTEGRATION_EPS_REL,
)

t0 = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY


def func_Ne(
    E: float, ka: float, q_Z: complex, m: int, n: int, s1: int, s2: int
) -> Callable[[float, float], float]:
    if s1 == s2:
        return lambda ka, q_Z: E
    elif s1 == 1 and s2 == 2:
        return lambda ka, q_Z: t0 + 2 * t0 * cos(q_Z) * exp(+1j * ka)
    elif s1 == 2 and s2 == 1:
        exponent_sign = 1 if (m == 0 and n == 0) else -1
        return lambda ka, q_Z: t0 + 2 * t0 * cos(q_Z) * exp(exponent_sign * 1j * ka)
    else:
        raise ValueError("Error! Check inputs s1 and s2 are valid. Must equal 1 or 2. ")


def green_function(Energy: float, m: int, n: int, s1: int, s2: int):
    def ka_integrand(ka, E, m, n, s1, s2):
        sum = 0

        for pm in [-1, 1]:
            q_Z = pm * acos(
                -0.5 * (cos(ka) + pm * sqrt((E**2 / t0**2) - sin(ka) * sin(ka)))
            )

            if q_Z.imag < 0:
                q_Z = -q_Z

            Ne = func_Ne(E, ka, q_Z, m, n, s1, s2)

            sum += (
                Ne(ka, q_Z)
                * exp(1j * (ka * (m + n) + q_Z * (m - n)))
                / (sin(2 * q_Z) + sin(q_Z) * cos(ka))
            )

        return sum

    Energy = Energy + 1j * ETA

    GF, _ = quad(
        ka_integrand,
        a=-pi,
        b=pi,
        args=(Energy, m, n, s1, s2),
        complex_func=True,
        limit=INTEGRATION_LIMIT,
        epsabs=INTEGRATION_EPS_ABS,
        epsrel=INTEGRATION_EPS_REL,
    )

    return 0.5 * (1j / (4.0 * pi * t0 * t0)) * GF
