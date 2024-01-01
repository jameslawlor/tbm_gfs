from scipy.integrate import nquad
from cmath import pi
from tbm_gfs.green_functions.graphene.functions import (
    Ne_lambda_function,
    dispersion_relation_lambda_function,
    phase_lambda_function,
)
from tbm_gfs.constants import (
    FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY,
    ETA,
)

from tbm_gfs.green_functions.graphene.config import (
    integration_options,
    kz_integration_bounds,
    ka_integration_bounds,
)
import numpy as np

t0 = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY


def green_function(Energy: float, m: int, n: int, s1: int, s2: int) -> complex:
    """
    Paul's thesis eq 3.1.13
    """

    def integrand(kz, ka, E, m, n, s1, s2):
        Ne = Ne_lambda_function(E, kz, ka, m, n, s1, s2)
        dispersion_epsilon = dispersion_relation_lambda_function(kz, ka)
        phase_term = phase_lambda_function(kz, ka, m, n)

        w = (Ne(kz, ka) * phase_term(kz, ka, m, n)) / (
            E**2 - dispersion_epsilon(kz, ka) ** 2
        )
        return w

    def real_integrand(kz, ka, E, m, n, s1, s2):
        return np.real(integrand(kz, ka, E, m, n, s1, s2))

    def imag_integrand(kz, ka, E, m, n, s1, s2):
        return np.imag(integrand(kz, ka, E, m, n, s1, s2))

    Energy = Energy + 1j * ETA

    GF_real, _ = nquad(
        real_integrand,
        [kz_integration_bounds, ka_integration_bounds],
        args=(Energy, m, n, s1, s2),
        opts=[integration_options, integration_options],
    )

    GF_imag, _ = nquad(
        imag_integrand,
        [kz_integration_bounds, ka_integration_bounds],
        args=(Energy, m, n, s1, s2),
        opts=[integration_options, integration_options],
    )

    GF = GF_real + 1j * GF_imag
    return (1 / (2 * pi * pi)) * GF
