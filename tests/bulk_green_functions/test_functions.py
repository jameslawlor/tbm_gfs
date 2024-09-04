import pytest

from tbm_gfs.bulk_green_functions.functions import (
    graphene_sector_method,
)

@pytest.mark.parametrize("m, n, s1, s2, expected", [
    (2, 1, 1, 1, (2, 1, 1, 1)),
    (5, 1, 1, 2, (5, 1, 1, 2)),
    (5, 5, 2, 2, (5, 5, 1, 1)),
    (0, 0, 2, 1, (0, 0, 1, 2)), # w-b case
])
def test_sector_0(m, n, s1, s2, expected):
    assert graphene_sector_method(m, n, s1, s2) == expected


@pytest.mark.parametrize("m, n, s1, s2, expected", [
    (1, 1, 1, 2, (1, 1, 2, 1)),
]) 
def test_sector_1(m, n, s1, s2, expected):
    assert graphene_sector_method(m, n, s1, s2) == expected

# def test_sector_2():
#     # Test for sector 2 where abs(n) > m
#     assert graphene_sector_method(2, -3, 0, 1) == (2, 0, 0, 1)

# def test_sector_3():
#     # Test for sector 3 where m < 0 and n < 0
#     assert graphene_sector_method(-2, -3, 0, 1) == (3, 2, 1, 0)

# def test_sector_4():
#     # Test for sector 4 where abs(m) > n
#     assert graphene_sector_method(-4, 1, 1, 0) == (1, 0, 1, 0)

# def test_sector_5():
#     # Test for sector 5 where abs(m) <= n
#     assert graphene_sector_method(-1, 4, 1, 0) == (1, 0, 0, 1)

# def test_irreducible_sector():
#     # Test for ensuring m >= n after transformation
#     assert graphene_sector_method(1, 4, 0, 1) == (4, 1, 0, 1)

# if __name__ == "__main__":
#     pytest.main()
