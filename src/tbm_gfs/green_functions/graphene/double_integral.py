from typing import Callable
from scipy.integrate import nquad
from cmath import cos, pi, sqrt, exp
from tbm_gfs.constants import (
    FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY,
    ETA,
    INTEGRATION_LIMIT,
    INTEGRATION_EPS_ABS,
    INTEGRATION_EPS_REL,
)
import numpy as np

t0 = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY


def func_Ne(
    E: float, kz: float, ka: float, m: int, n: int, s1: int, s2: int
) -> Callable[[float, float], float]:
    if s1 == s2:
        return lambda kz, ka: E
    elif s1 == 1 and s2 == 2:
        return lambda kz, ka: t0 * (1 + 2 * cos(kz) * exp(-1j * ka))
    elif s1 == 2 and s2 == 1:
        exponent_sign = -1 if (m == 0 and n == 0) else 1
        return lambda kz, ka: t0 * (1 + 2 * cos(kz) * exp(exponent_sign * 1j * ka))
    else:
        raise ValueError("Error! Check inputs s1 and s2 are valid. Must equal 1 or 2. ")


def green_function(Energy: float, m: int, n: int, s1: int, s2: int) -> complex:
    """
    Paul's thesis 3.1.13
    """

    def integrand(kz, ka, E):
        dispersion_epsilon = t0 * sqrt(
            1 + 4 * cos(ka) * cos(kz) + 4 * cos(kz) * cos(kz)
        )

        Ne = func_Ne(E, kz, ka, m, n, s1, s2)

        return (Ne(kz, ka) * exp(1j * (ka * (m + n) + kz * (m - n)))) / (
            E**2 - dispersion_epsilon**2
        )

    def real_integrand(kz, ka, E):
        return np.real(integrand(kz, ka, E))

    def imag_integrand(kz, ka, E):
        return np.imag(integrand(kz, ka, E))

    Energy = Energy + 1j * ETA

    options = {
        "limit": INTEGRATION_LIMIT,
        "epsabs": INTEGRATION_EPS_ABS,
        "epsrel": INTEGRATION_EPS_REL,
    }

    GF_real, _ = nquad(
        real_integrand,
        [[-pi / 2.0, pi / 2.0], [-pi, pi]],
        args=(Energy,),
        opts=[options, options],
    )

    GF_imag, _ = nquad(
        imag_integrand,
        [[-pi / 2.0, pi / 2.0], [-pi, pi]],
        args=(Energy,),
        opts=[options, options],
    )

    GF = GF_real + 1j * GF_imag
    return (1 / (2 * pi * pi)) * GF
