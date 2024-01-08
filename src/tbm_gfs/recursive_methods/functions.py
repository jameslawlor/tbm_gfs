import numpy as np

from tbm_gfs.constants import (
    ETA,
)

def dyson_eq(
    cell_to_connect_gf,
    current_surface_gf,
    cell_to_surface_connection,
    surface_to_cell_connection,
    identity_matrix,
):
    """
    See Section 2.5 in my thesis
    """
    return np.dot(
                np.linalg.inv(
                    identity_matrix
                    - np.dot(
                        np.dot(cell_to_connect_gf, cell_to_surface_connection),
                        np.dot(current_surface_gf, surface_to_cell_connection ),
                    )
                ),
                cell_to_connect_gf,
            )

def check_convergence(matrix_i, matrix_i_minus_1) -> bool:
    """
    Checks if matrix_i and matrix_i_minus_1
    are sufficiently similar we consider we have
    achieved convergence - the matrices then
    represent surface unit cells in a semi-infinite
    lead structure
    """
    return np.allclose(
        matrix_i,
        matrix_i_minus_1,
        rtol=1e-04,
        atol=1e-04,
    )


def standard_recursive_lead_generator(
    unit_cell_hamiltonian, energy, connection_matrix_LR, connection_matrix_RL
):
    # TODO: Generalise function to work with opposite direction lead 
    
    identity_matrix = np.eye(unit_cell_hamiltonian.shape[0], dtype=complex)
    green_function_unit_cell = np.linalg.inv(
        (energy + 1j * ETA) * identity_matrix - unit_cell_hamiltonian
    )

    surface_lead_gf = green_function_unit_cell  # initialise
    surface_lead_gf_previous_iteration = None
    is_converged= False

    while not is_converged:

        surface_lead_gf_previous_iteration = surface_lead_gf

        surface_lead_gf = dyson_eq(
            cell_to_connect_gf = green_function_unit_cell,
            current_surface_gf = surface_lead_gf_previous_iteration,
            cell_to_surface_connection = connection_matrix_LR,
            surface_to_cell_connection = connection_matrix_RL,
            identity_matrix = identity_matrix,
        )
        is_converged = check_convergence(surface_lead_gf, surface_lead_gf_previous_iteration)

    return surface_lead_gf


if __name__ == "__main__":
    from tbm_gfs.recursive_methods.matrices.zigzag_nanoribbon import (
        unit_cell_hamiltonian,
        LR_connection_matrix,
        RL_connection_matrix,
    )
    ldos_list = []
    e_range = np.linspace(-4,4,1001)
    for e in e_range:
        size = 4
        H_UC = unit_cell_hamiltonian(size)
        V_LR = LR_connection_matrix(size)
        V_RL = RL_connection_matrix(size)

        S_L = standard_recursive_lead_generator(H_UC,e,V_LR,V_RL )
        S_R = standard_recursive_lead_generator(H_UC,e,V_RL,V_LR )

        ident = np.eye(S_L.shape[0], dtype=complex)
        G_final = dyson_eq(
            cell_to_connect_gf=S_L,
            current_surface_gf=S_R,
            cell_to_surface_connection=V_LR,
            surface_to_cell_connection=V_RL,
            identity_matrix=ident,
        )
        ldos = G_final[3][3].imag
        ldos_list.append(ldos)
        print(e,ldos)

    import matplotlib.pylab as plt
    plt.plot(e_range, ldos_list)
    plt.show()


# def left_lead_recursive(Energy, H00):
# 	g00 = invert((Energy+i*eta)*I - H00) 
# 	SL = g00
# 	SLold = np.zeros(shape=(2*size+1,2*size+1))
# 	#Build Left Lead
# 	while np.abs(np.sum(SL - SLold)) > tolerance:
# 		SLold = np.array(SL)
# 		SL = np.dot( invert(I - np.dot(np.dot(g00,VRL),np.dot(SLold,VLR) )) , g00)
# 	return SL