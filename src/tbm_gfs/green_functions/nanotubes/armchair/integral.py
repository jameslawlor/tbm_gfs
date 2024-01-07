from scipy.integrate import quad
from cmath import pi

from tbm_gfs.green_functions.functions import (
    Ne_lambda_function,
    phase_lambda_function,
    dispersion_relation_lambda_function,
)

from tbm_gfs.constants import (
    INTEGRATION_LIMIT,
    INTEGRATION_EPS_ABS,
    INTEGRATION_EPS_REL,
    ETA,
    FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY,
)

integration_options = {
    "limit": INTEGRATION_LIMIT,
    "epsabs": INTEGRATION_EPS_ABS,
    "epsrel": INTEGRATION_EPS_REL,
}

t0 = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY


def integrand(kz, E, m, n, s1, s2, j, n_c):
    ka = pi * j / n_c
    Ne = Ne_lambda_function(E, kz, ka, m, n, s1, s2)
    phase_term = phase_lambda_function(kz, ka, m, n)
    dispersion_relation = dispersion_relation_lambda_function(kz, ka)

    return (
        Ne(kz, ka)
        * phase_term(kz, ka, m, n)
        / (E * E - dispersion_relation(kz, ka) ** 2)
    )


def green_function(
    n_c: int,
    Energy: float,
    m: int,
    n: int,
    s1: int,
    s2: int,
):
    Energy = Energy + 1j * ETA

    GF = 0
    for j in range(0, n_c):
        res, _ = quad(
            integrand,
            a=-pi,
            b=pi,
            args=(Energy, m, n, s1, s2, j, n_c),
            complex_func=True,
            limit=integration_options["limit"],
            epsabs=integration_options["epsabs"],
            epsrel=integration_options["epsrel"],
        )
        GF += res

    return GF / (2 * pi * n_c)
