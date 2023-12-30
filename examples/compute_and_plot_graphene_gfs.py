from tbm_gfs.green_functions.graphene.single_integral import green_function
import numpy as np
import matplotlib.pylab as plt

m = 1
n = 1
s1 = 1
s2 = 2

energy_range = np.arange(-3.5, 3.5, 0.025)
gfs = [green_function(energy, m, n, s1, s2) for energy in energy_range]
plt.plot(
    energy_range,
    np.real(gfs),
    color="green",
)
plt.plot(
    energy_range,
    np.imag(gfs),
    color="red",
)

m = -1
n = -1
s1 = 2
s2 = 1
gfs = [green_function(energy, m, n, s1, s2) for energy in energy_range]
plt.plot(
    energy_range,
    np.real(gfs),
    color="green",
    linestyle='dashed',
    marker='o'
)
plt.plot(
    energy_range,
    np.imag(gfs),
    color="red",
    linestyle='dashed',
    marker='o'
)

plt.show()
