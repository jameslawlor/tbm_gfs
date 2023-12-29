from scipy.integrate import quad
from cmath import sin, cos, acos, pi


def g00(*args):
    t0 = -1.0

    def kz_integrand(kz, E):
        if E.imag == 0.0:
            E += 1j * 1.0e-3

        # the complex pole
        q_A = acos(
            (E**2.0 - (t0 * t0) - 4.0 * t0 * t0 * (cos(kz) ** 2.0))
            / (4.0 * t0 * t0 * cos(kz))
        )

        if q_A.imag < 0:
            q_A = -q_A

        W = E / (cos(kz) * sin(q_A))
        return W

    GF, _ = quad(
        kz_integrand,
        a=-pi / 2.0,
        b=pi / 2.0,
        args=args,
        complex_func=True,
        limit=300
    )
    return (1j / (4.0 * pi * t0 * t0)) * GF
