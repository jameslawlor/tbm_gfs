from tbm_gfs.green_functions.graphene.single_integral import green_function
import numpy as np
import matplotlib.pylab as plt

energy_range = np.arange(-4.0, 4.0, 0.2)
gfs = [green_function(energy) for energy in energy_range]
plt.plot(energy_range, np.imag(gfs))
plt.show()
