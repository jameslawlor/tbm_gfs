"""
See Equation 13 in
Duffy, J. M., et al.
"Variable range of the RKKY interaction in edged graphene."
Journal of Physics: Condensed Matter 26.5 (2013): 055007.
"""

from scipy.integrate import quad
from cmath import sin, cos, acos, pi, exp
from typing import Callable
import numpy as np

from tbm_gfs.constants import (
    FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY,
    ETA,
    INTEGRATION_LIMIT,
    INTEGRATION_EPS_ABS,
    INTEGRATION_EPS_REL,
)

t0 = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY


def func_Ne(
    E: float, kz: float, q_A: complex, m: int, n: int, s1: int, s2: int
) -> Callable[[float, float], float]:
    if s1 == s2:
        return lambda kz, q_A: E
    elif s1 == 1 and s2 == 2:
        return lambda kz, q_A: t0 + 2 * t0 * cos(kz) * exp(+1j * q_A)
    elif s1 == 2 and s2 == 1:
        exponent_sign = 1 if (m == 0 and n == 0) else -1
        return lambda kz, q_A: t0 + 2 * t0 * cos(kz) * exp(exponent_sign * 1j * q_A)
    else:
        raise ValueError("Error! Check inputs s1 and s2 are valid. Must equal 1 or 2. ")


def green_function(Energy: float, m: int, n: int, s1: int, s2: int):
    def kz_integrand(kz, E, m, n, s1, s2):
        # the complex pole
        q_A = acos(
            (E**2.0 - (t0 * t0) - 4.0 * t0 * t0 * (cos(kz) ** 2.0))
            / (4.0 * t0 * t0 * cos(kz))
        )

        if q_A.imag < 0:
            q_A = -q_A

        Ne = func_Ne(E, kz, q_A, m, n, s1, s2)

        return (
            Ne(kz, q_A)
            * exp(1j * (q_A * (m + n) + kz * (m - n)))
            / (cos(kz) * sin(q_A))
        )

    Energy = Energy + 1j * ETA
    
    GF, _ = quad(
        kz_integrand,
        a=-pi / 2.0,
        b=pi / 2.0,
        args=(Energy, m, n, s1, s2),
        complex_func=True,
        limit=INTEGRATION_LIMIT,
        epsabs=INTEGRATION_EPS_ABS,
        epsrel=INTEGRATION_EPS_REL,
    )

    return (1j / (4.0 * pi * t0 * t0)) * GF