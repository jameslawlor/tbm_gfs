from tbm_gfs.bulk_green_functions.graphene.single_integral import green_function
from cmath import pi


def local_density_of_states(
    lam: float, Energy: float, m: int, n: int, s1: int, s2: int
):
    """
    Compute variations in local density of states
    in graphene at a specified lattice location `i`
    caused by a perturbation at `A` represented
    by an energy `lam`.
    """

    gia = green_function(Energy=Energy, m=m, n=n, s1=s1, s2=s2)
    gaa = green_function(Energy=Energy, m=0, n=0, s1=1, s2=1)

    dG = gia * gia * lam / (1 - gaa * lam)

    return (-1.0 / pi) * dG.imag


if __name__ == "__main__":
    from tbm_gfs.plotter import GreenFunctionPlotter

    energy = 0.2
    s1 = 1
    s2 = 2
    distances = range(10, 80)
    lam = 0.1

    delta_ldos_same_sublattice = [
        local_density_of_states(lam, energy, m, m, s1, s1) for m in distances
    ]
    delta_ldos_opp_sublattice = [
        local_density_of_states(lam, energy, m, m, s1, s2) for m in distances
    ]

    import matplotlib.pylab as plt

    plt.figure()
    plt.xlabel("m")
    plt.ylabel("Change in LDOS")
    plt.plot(
        distances,
        delta_ldos_same_sublattice,
        label="same sublattice",
        color="k",
        marker="o",
    )
    plt.plot(
        distances,
        delta_ldos_opp_sublattice,
        label="opp sublattice",
        color="r",
        marker="o",
    )
    plt.legend()
    plt.grid(True)
    plt.show()
