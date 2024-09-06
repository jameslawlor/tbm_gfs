import pytest

from tbm_gfs.bulk_green_functions.functions import (
    graphene_sector_method,
)


@pytest.mark.parametrize(
    "m, n, s1, s2, expected",
    [
        (2, 1, 1, 1, (2, 1, 1, 1)),
        (5, 1, 1, 2, (5, 1, 1, 2)),
        (5, 5, 2, 2, (5, 5, 1, 1)),
        (0, 0, 2, 1, (0, 0, 1, 2)),  # w-b case
    ],
)
def test_sector_0(m, n, s1, s2, expected):
    assert graphene_sector_method(m, n, s1, s2) == expected


@pytest.mark.parametrize(
    "m, n, s1, s2, expected",
    [
        (1, -1, 1, 2, (1, 1, 2, 1)),
        (2, -1, 1, 1, (1, 1, 1, 1)),
        (3, -1, 2, 1, (1, 1, 1, 2)),
        (3, -2, 1, 1, (2, 1, 1, 1)),
    ],
)
def test_sector_1(m, n, s1, s2, expected):
    result = graphene_sector_method(m, n, s1, s2)
    assert result == expected


@pytest.mark.parametrize(
    "m, n, s1, s2, expected",
    [
        (1, -2, 1, 1, (1, 1, 1, 1)),
        (1, -2, 1, 2, (1, 0, 1, 2)),
        (1, -2, 2, 1, (2, 1, 2, 1)),
        (1, -3, 1, 1, (2, 1, 1, 1)),
    ],
)
def test_sector_2(m, n, s1, s2, expected):
    result = graphene_sector_method(m, n, s1, s2)
    assert result == expected


@pytest.mark.parametrize(
    "m, n, s1, s2, expected",
    [
        (-1, 0, 1, 1, (1, 0, 1, 1)),
        (-1, -1, 1, 1, (1, 1, 1, 1)),
        (-1, -1, 1, 2, (1, 1, 2, 1)),
        (-2, -1, 1, 1, (2, 1, 1, 1)),
    ],
)
def test_sector_3(m, n, s1, s2, expected):
    result = graphene_sector_method(m, n, s1, s2)
    assert result == expected


@pytest.mark.parametrize(
    "m, n, s1, s2, expected",
    [
        (-1, 1, 1, 1, (1, 0, 1, 1)),
        (-2, 1, 1, 1, (1, 1, 1, 1)),
        (-3, 1, 1, 1, (2, 1, 1, 1)),
        (-2, 0, 1, 2, (2, 0, 2, 1)),
        (-1, 0, 2, 1, (1, 0, 1, 2)),
    ],
)
def test_sector_4(m, n, s1, s2, expected):
    result = graphene_sector_method(m, n, s1, s2)
    assert result == expected


@pytest.mark.parametrize(
    "m, n, s1, s2, expected",
    [
        (-1, 2, 1, 1, (1, 1, 1, 1)),
        (-1, 3, 1, 1, (2, 1, 1, 1)),
        (-1, 1, 1, 2, (1, 1, 2, 1)),
        (-1, 1, 2, 1, (1, 1, 2, 1)),
    ],
)
def test_sector_5(m, n, s1, s2, expected):
    result = graphene_sector_method(m, n, s1, s2)
    assert result == expected
