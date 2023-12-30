import pytest
import numpy as np
from tbm_gfs.green_functions.graphene.single_integral import (
    green_function as single_integral_gf,
)

MAX_REL_ERROR = 2e-2 # Allow only 2% error comparing double vs single integrals

@pytest.fixture
def sample_energies():
    incr = 0.25
    e_max = 4.0
    return np.arange(-e_max, e_max, incr)

def test_g00_bw_wb(sample_energies):
    result_1 = single_integral_gf(sample_energies, m=0, n=0, s1=1, s2=2)
    result_2 = single_integral_gf(sample_energies, m=0, n=0, s1=2, s2=1)
    assert result_1.real ==  pytest.approx(
            result_2.real, rel=MAX_REL_ERROR
        )
    assert result_1.imag ==  pytest.approx(
            result_2.imag, rel=MAX_REL_ERROR
        )