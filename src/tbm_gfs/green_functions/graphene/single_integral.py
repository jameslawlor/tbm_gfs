from scipy.integrate import quad
from cmath import sin, cos, acos, pi
from tbm_gfs.constants import (
    FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY,
    ETA,
    INTEGRATION_LIMIT,
)

t0 = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY


def green_function(Energy):
    """
    See Equation 13 in
    Duffy, J. M., et al.
    "Variable range of the RKKY interaction in edged graphene."
    Journal of Physics: Condensed Matter 26.5 (2013): 055007.
    """

    def kz_integrand(kz, E):
        # the complex pole
        q_A = acos(
            (E**2.0 - (t0 * t0) - 4.0 * t0 * t0 * (cos(kz) ** 2.0))
            / (4.0 * t0 * t0 * cos(kz))
        )

        if q_A.imag < 0:
            q_A = -q_A

        return E / (cos(kz) * sin(q_A))

    if Energy.imag == 0.0:
        Energy += 1j * ETA

    GF, _ = quad(
        kz_integrand,
        a=-pi / 2.0,
        b=pi / 2.0,
        args=(Energy),
        complex_func=True,
        limit=INTEGRATION_LIMIT,
    )
    return (1j / (4.0 * pi * t0 * t0)) * GF
