"""
Stationary Phase Approximation (SPA) is a closed-form approximation of the graphene Green Functions (single or double integral) solutions
Works best in high-symmetry lattice directions (armchair / zigzag) over larger distances (i.e. m+n >> 0) and at Energies close to the Dirac point
Basic idea is to use a Taylor expansion (every physicist's best friend) 
Originally derived by Stephen, but my derivation is in Section 3.2.2 of my thesis "Electronic properties of doped carbon-based nanostructures"
Google Scholar link: https://scholar.google.com/citations?view_op=view_citation&hl=en&user=XeiCaFIAAAAJ&authuser=1&citation_for_view=XeiCaFIAAAAJ:geHnlv5EZngC
""" 

from cmath import sin, cos, acos, pi, sqrt, exp

from tbm_gfs.constants import (
    FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY,
    ETA,
)

from tbm_gfs.green_functions.graphene.functions import (
    Ne_lambda_function,
    phase_lambda_function,
)

t0 = FIRST_NEAREST_NEIGHBOUR_HOPPING_ENERGY

def spa_same_sublattice_armchair_gf(E,m,n,s1,s2):

	D=m+n
	# # # Eq 3.16 in my thesis
	if abs(E) < abs(t0):
		kz=0
		qa = (-E/abs(E))*acos(-sqrt(1-(E**2 / t0**2)))
		Ne = Ne_lambda_function(E, kz, qa, m, n, s1, s2)
		phase_term = phase_lambda_function(kz, qa, m, n)
		denom = sqrt(E*(E**2 + 3*t0**2)*sqrt(t0**2 - E**2))
		return - Ne(kz,qa)*sqrt((2*1j)/(pi*D))*phase_term(kz,qa,m,n)/denom
	else:

		from math import copysign
		sig = copysign(1.0,-E.real)

		qa = sig*acos((E**2 - 5*t0**2) / (4*t0**2))
		kz = -acos((sqrt(t0**2 - E**2)/(2*t0)))
		Ne = Ne_lambda_function(E, kz, qa, m, n, s1, s2)
		phase_term = lambda kz, ka, m, n: exp(1j * (ka * (m + n)))
		denom = sqrt(E**2 + 3*t0**2)*(((t0**2 - E**2)*(E**2 - 9*t0**2))**(1/4))
		
		return sig*1j*Ne(kz,qa)*sqrt(-sig*(2*1j)/(pi*D))*phase_term(kz,qa,m,n)/denom