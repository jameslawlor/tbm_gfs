from tbm_gfs.bulk_green_functions.graphene.single_integral import green_function
from tbm_gfs.bulk_green_functions.graphene.config import integration_options
from tbm_gfs.constants import ETA

from cmath import pi
from scipy.integrate import quad

def charge_density(lam:float, fermi_energy: float, m: int, n: int, s1: int, s2: int):
    """
    Compute variations in charge density
    in graphene at a specified lattice location `i`
    caused by a perturbation at `A` represented
    by an energy `lam`. 
    """
    def _integrand(s, lam , Ef, m , n , s1 , s2):
        Energy = Ef + 1j*((1.0 + ETA - s)/s)
        gia = green_function(Energy=Energy, m=m, n=n, s1=s1, s2=s2)
        gaa = green_function(Energy=Energy, m=0, n=0, s1=1, s2=1)
        return (1/(s*s)) * gia*gia*lam / (1.0 - gaa*lam)

    dn, _ = quad(
            _integrand,
            a=0,
            b=1,
            args=(lam, fermi_energy, m, n, s1, s2),
            complex_func=True,
            limit=integration_options["limit"],
            epsabs=integration_options["epsabs"],
            epsrel=integration_options["epsrel"],
        )

    return (2.0/pi) * dn.real


if __name__ == '__main__':
    from tbm_gfs.plotter import GreenFunctionPlotter

    fermi_energy = 0.2
    s1 = 1
    s2 = 2
    distances = range(10,40)
    lam = -1
    
    delta_charge_density_same_sublattice_sites = [charge_density(lam, fermi_energy, m, m, s1, s1) for m in distances]
    delta_charge_density_opp_sublattice_sites = [charge_density(lam, fermi_energy, m, m, s1, s2) for m in distances]
    import matplotlib.pylab as plt

    plt.figure()
    plt.xlabel("m")
    plt.ylabel("Change in charge density")
    plt.plot(distances, delta_charge_density_same_sublattice_sites, label='same sublattice', color='k', marker='o')
    plt.plot(distances, delta_charge_density_opp_sublattice_sites, label='opp sublattice', color='r', marker='o')
    plt.legend()
    plt.grid(True)
    plt.show()