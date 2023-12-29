from scipy.integrate import dblquad
from cmath import cos, pi, sqrt
from tbm_gfs.constants import FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY, ETA
import numpy as np

t0 = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY


def green_function(Energy):
    """
    Paul's thesis 3.1.13
    """

    def integrand(kz, ka, E):
        dispersion_epsilon = t0 * sqrt(
            1 + 4 * cos(ka) * cos(kz) + 4 * cos(kz) * cos(kz)
        )

        return E / (E**2 - dispersion_epsilon**2)

    def real_integrand(kz, ka, E):
        return np.real(integrand(kz, ka, E))

    def imag_integrand(kz, ka, E):
        return np.imag(integrand(kz, ka, E))

    if Energy.imag == 0.0:
        Energy += 1j * ETA

    GF_real, _ = dblquad(
        real_integrand,
        a=-pi / 2.0,
        b=pi / 2.0,
        gfun=-pi,
        hfun=pi,
        args=(Energy,),
        epsabs=1e-5,
        epsrel=1e-5,
    )

    GF_imag, _ = dblquad(
        imag_integrand,
        a=-pi / 2.0,
        b=pi / 2.0,
        gfun=-pi,
        hfun=pi,
        args=(Energy,),
        epsabs=1e-5,
        epsrel=1e-5,
    )

    GF = GF_real + 1j * GF_imag
    return (1 / (2 * pi * pi)) * GF
