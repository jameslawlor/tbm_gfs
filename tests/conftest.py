import pytest


@pytest.fixture
def max_relative_integration_comparison_error():
    return 2e-2  # Allow only 2% difference when comparing test results of two numerical integrals
