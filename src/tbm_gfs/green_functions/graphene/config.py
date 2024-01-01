from tbm_gfs.constants import (
    INTEGRATION_LIMIT,
    INTEGRATION_EPS_ABS,
    INTEGRATION_EPS_REL,
)
from cmath import pi

integration_options = {
    "limit": INTEGRATION_LIMIT,
    "epsabs": INTEGRATION_EPS_ABS,
    "epsrel": INTEGRATION_EPS_REL,
}

kz_integration_bounds = (-pi / 2.0, pi / 2.0)
ka_integration_bounds = (-pi, pi)
