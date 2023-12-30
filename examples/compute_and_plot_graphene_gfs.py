from tbm_gfs.green_functions.graphene.single_integral import green_function
import numpy as np
import matplotlib.pylab as plt

m = 40
n = 40
s1 = 1
s2 = 1

energy_range = np.arange(-0.5, 0.5, 0.01)
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
plt.show()
