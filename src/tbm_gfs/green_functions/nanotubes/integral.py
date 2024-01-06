from scipy.integrate import quad
from cmath import pi, exp, cos, sin, acos, sqrt

from tbm_gfs.constants import (
    INTEGRATION_LIMIT,
    INTEGRATION_EPS_ABS,
    INTEGRATION_EPS_REL,
    ETA,
)
import numpy as np
import matplotlib.pylab as plt

integration_options = {
    "limit": INTEGRATION_LIMIT,
    "epsabs": INTEGRATION_EPS_ABS,
    "epsrel": INTEGRATION_EPS_REL,
}


def integrand(kz, E, m, n, s1, s2, j, n_c):
    f = 1 + 2 * cos(kz) * exp(1j * pi * j / n_c)
    return (
        E
        * exp(1j * (pi * j / n_c) * (m + n) + 1j * kz * (m - n))
        / (E * E - f * np.conj(f))
    )


def armchair_nanotube_green_function(
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


t = -1
i = 1j


def acnt(En, Nc, m, n, S1, S2):
    En+=0.001*i
    if (m + n) < 0:
        S1, S2 = S2, S1
        m, n = -m, -n
    G = 0
    for j in range(1, Nc + 1):
        kA = pi * j / Nc
        for w in (-1, 1):
            q = acos(
                (-0.5) * (t / t) * cos(kA)
                + w * (0.5 / t) * sqrt(En * En - (t**2) * sin(kA) ** 2)
            )
            if q.imag < 0:
                q = -q

            if S1 == S2:
                term = En
            elif S1 < S2:
                term = t + 2 * t * exp(i * kA) * cos(q)
            else:
                term = t + 2 * t * exp(-i * kA) * cos(q)
            G += (
                term
                * (
                    ((i) / (4 * Nc * t))
                    * exp(i * kA * (m + n))
                    * exp(i * (q * abs(m - n)))
                )
                / (sin(q) * (t * cos(kA) + 2 * t * cos(q)))
            )
    return G


def zznt(En, Nc, m, n, S1, S2):
    if (m + n) < 0:
        S1, S2 = S2, S1
        m, n = -m, -n
    G = 0
    for j in range(1, Nc + 1):
        kZ = j * pi / Nc
        q = acos((En**2 - t**2 - 4 * (t**2) * (cos(kZ) ** 2)) / (4 * t * t * cos(kZ)))
        if q.imag < 0:
            q = -q

        if (j == Nc / 2.0) and (m + n == 0):
            if S1 == S2:
                Ne = En
            else:
                Ne = t
            G += (Ne * exp(i * kZ * (m - n))) / (Nc * (En**2 - t**2))
        else:
            if S1 == S2:
                Ne = En
            if S1 < S2:
                Ne = t + 2 * t * cos(kZ) * exp(+i * q)
            if S1 > S2:
                if m == 0 and n == 0:
                    Ne = t + 2 * t * cos(kZ) * exp(+i * q)
                else:
                    Ne = t + 2 * t * cos(kZ) * exp(-i * q)
            G += (i / (4 * Nc * (t * t))) * (
                (Ne * exp(i * 2 * kZ * (m - n)) * exp(i * q * (m + n)))
                / (cos(kZ) * sin(q))
            )
    return G


if __name__ == "__main__":
    en = np.arange(-3, 3, 0.01)
    gfs = [armchair_nanotube_green_function(11, e, 0, 0, 1, 1) for e in en]
    gfs2 = [acnt(e, 11, 0, 0, 1, 1) for e in en]
    plt.plot(
        en,
        [gf.real for gf in gfs],
        color="red",
    )
    plt.plot(en, [gf.real for gf in gfs2], color="black", marker="o", linestyle='None')
    plt.show()
    plt.plot(
        en,
        [gf.imag for gf in gfs],
        color="red",
    )
    plt.plot(en, [gf.imag for gf in gfs2], color="black", marker="o", linestyle='None')
    plt.show()
