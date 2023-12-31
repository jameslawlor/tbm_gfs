# Typically -2.7eV. See Reich, Stephanie, et al. "Tight-binding description of graphene." Physical Review B 66.3 (2002): 035412.
FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY = -1.0
ETA = 1.0e-3  # small imaginary part added to energies during integration, to prevents poles causing divergencies
INTEGRATION_LIMIT = 100  # for scipy.integrate.quad
INTEGRATION_EPS_ABS = 1e-4
INTEGRATION_EPS_REL = 1e-4
