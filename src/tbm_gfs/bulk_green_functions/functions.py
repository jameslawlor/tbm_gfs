from typing import Callable
from cmath import cos, exp, sqrt

from tbm_gfs.constants import (
    FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY,
)

t0 = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY

# TODO: Refactor lambda functions to clean up unneeded arguments


def Ne_lambda_function(
    E: complex, kz: complex, ka: complex, m: int, n: int, s1: int, s2: int
) -> Callable[[complex, complex], complex]:
    if s1 == s2:
        return lambda kz, ka: E
    elif s1 == 1 and s2 == 2:
        return lambda kz, ka: t0 * (1 + 2 * cos(kz) * exp(1j * ka))
    elif s1 == 2 and s2 == 1:
        exponent_sign = 1 if (m == 0 and n == 0) else -1
        return lambda kz, ka: t0 * (1 + 2 * cos(kz) * exp(exponent_sign * 1j * ka))
    else:
        raise ValueError("Error! Check inputs s1 and s2 are valid. Must equal 1 or 2. ")


def phase_lambda_function(
    kz: complex,
    ka: complex,
    m: int,
    n: int,
) -> Callable[[complex, complex, int, int], complex]:
    return lambda kz, ka, m, n: exp(1j * (ka * (m + n) + kz * (m - n)))


def dispersion_relation_lambda_function(
    kz: complex, ka: complex
) -> Callable[[complex, complex], complex]:
    return lambda kz, ka: t0 * sqrt(1 + 4 * cos(ka) * cos(kz) + 4 * cos(kz) * cos(kz))

    
def graphene_sector_method(m, n, s1, s2):
    """
	Transform any (m,n,s1,s2) vector on the graphene lattice
	into an equivalent one in the irreducible sector from angle 0 to pi/6
    
    See Paul's thesis Appendix A
    """

    s = s2 - s1

    # Handle white-black case
    if s < 0:
        m, n, s = -n, -m, -s

    # Determine the correct sector and adjust m, n, s1, s2 accordingly
    if m >= 0 and n >= 0:
        pass  # Sector 0: no change needed
    elif m >= 0 and abs(n) <= m:  # Sector 1 - bottom right 
        m, n, s = -n, m + n + s, -s
    elif m >= 0 and abs(n) > m:  # Sector 2 - bottom left
        m, n, s = m, -n - m - s, s
    elif m < 0 and n < 0:  # Sector 3 - left
        m, n, s = -n, -m, -s
    elif abs(m) > n:  # Sector 4 - upper left
        m, n, s = n, -n -m -s, s
    elif m < 0 and abs(m) <= 0:  # Sector 5 - upper right
        m, n, s = -m, n + m + s, -s

    # Ensure m >= n to move to the irreducible sector
    if n > m:
        m, n = n, m

    if not (0 <= n and n <= m):
        raise ValueError("No")

    if s == 0:
        return m, n, 1, 1
    elif s > 0:
        return m, n, 1, 2
    elif s < 0:
        return m, n, 2, 1
