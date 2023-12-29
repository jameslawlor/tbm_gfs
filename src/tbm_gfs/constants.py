# Typically -2.7eV. See Reich, Stephanie, et al. "Tight-binding description of graphene." Physical Review B 66.3 (2002): 035412.
FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY = -1.0
ETA = 1.0e-4  # small imaginary part added to energies during integration, to prevents poles causing divergencies
INTEGRATION_LIMIT = 300  # for scipy.integrate.quad
