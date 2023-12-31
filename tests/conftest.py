import pytest


@pytest.fixture
def max_relative_integration_comparison_error():
    return 1e-3  # Allow only 0.1% difference when comparing test results of two numerical integrals
