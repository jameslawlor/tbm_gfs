from tbm_gfs.bulk_green_functions.graphene.single_integral import green_function
from cmath import pi

def ldos(lam:float, Energy: float, m: int, n: int, s1: int, s2: int):
    """
    Compute variations in local density of states
    in graphene at a specified lattice location `i`
    caused by a perturbation at `A` represented
    by an energy `lam`. 
    """

    gia = green_function(Energy=Energy, m=m, n=n, s1=s1, s2=s2)
    gaa = green_function(Energy=Energy, m=0, n=0, s1=1, s2=1)

    dG = gia * gia * lam / (1 - gaa * lam)

    return (-1.0/pi) * dG.imag


if __name__ == '__main__':
    import matplotlib.pylab as plt
    xrange=range(1,50)
    tmp = [ldos(-1, 0.2, m, m, 1, 1) for m in xrange]
    plt.plot(xrange,tmp,color='green')

    tmp = [ldos(-1, 0.2, m, m, 1, 2) for m in xrange]
    plt.plot(xrange,tmp,color='red')
    plt.show()