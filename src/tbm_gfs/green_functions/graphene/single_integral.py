"""
See Equation 13 in
Duffy, J. M., et al.
"Variable range of the RKKY interaction in edged graphene."
Journal of Physics: Condensed Matter 26.5 (2013): 055007.
"""

from scipy.integrate import quad
from cmath import sin, cos, acos, pi, sqrt
from tbm_gfs.green_functions.functions import (
    Ne_lambda_function,
    phase_lambda_function,
)
from tbm_gfs.green_functions.graphene.config import (
    integration_options,
    kz_integration_bounds,
    ka_integration_bounds,
)

from tbm_gfs.constants import (
    FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY,
    ETA,
)

t0 = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY


def ka_integrand(ka, E, m, n, s1, s2):
    sum = 0

    for pm in [-1, 1]:
        q_Z = pm * acos(
            -0.5 * (cos(ka) + pm * sqrt((E**2 / t0**2) - sin(ka) * sin(ka)))
        )

        if q_Z.imag < 0:
            q_Z = -q_Z

        Ne = Ne_lambda_function(E, q_Z, ka, m, n, s1, s2)
        phase_term = phase_lambda_function(q_Z, ka, m, n)

        sum += (
            Ne(q_Z, ka)
            * phase_term(q_Z, ka, m, n)
            / (sin(2 * q_Z) + sin(q_Z) * cos(ka))
        )

    return 0.5 * sum


def kz_integrand(kz, E, m, n, s1, s2):
    # the complex pole
    q_A = acos(
        (E**2.0 - (t0 * t0) - 4.0 * t0 * t0 * (cos(kz) ** 2.0))
        / (4.0 * t0 * t0 * cos(kz))
    )

    if q_A.imag < 0:
        q_A = -q_A

    Ne = Ne_lambda_function(E, kz, q_A, m, n, s1, s2)
    phase_term = phase_lambda_function(kz, q_A, m, n)

    return Ne(kz, q_A) * phase_term(kz, q_A, m, n) / (cos(kz) * sin(q_A))


def setup_integration(integration_variable):
    if integration_variable == "kz":
        return kz_integrand, kz_integration_bounds
    elif integration_variable == "ka":
        return ka_integrand, ka_integration_bounds
    else:
        raise ValueError(
            "Choice of integration_variable not recognised. Choose either ka or kz"
        )


def green_function(
    Energy: float, m: int, n: int, s1: int, s2: int, integration_variable="kz"
):
    Energy = Energy + 1j * ETA

    integrand, bounds = setup_integration(integration_variable)

    GF, _ = quad(
        integrand,
        a=bounds[0],
        b=bounds[1],
        args=(Energy, m, n, s1, s2),
        complex_func=True,
        limit=integration_options["limit"],
        epsabs=integration_options["epsabs"],
        epsrel=integration_options["epsrel"],
    )

    return GF * 1j / (4.0 * pi * t0 * t0)
