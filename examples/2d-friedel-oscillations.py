from tbm_gfs.bulk_green_functions.functions import (
    convert_lattice_vectors_to_real_space,
    graphene_sector_method,
)
from tbm_gfs.friedel_oscillations.graphene.ldos import local_density_of_states
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata


def draw_hexagonal_lattice(
    x_values,
    y_values,
    sublattice_values,
):
    # Define the nearest neighbors in a hexagonal lattice
    neighbor_offsets = [
        [1, 0],
        [-0.5, np.sqrt(3) / 2],
        [-0.5, -np.sqrt(3) / 2],
    ]

    # Draw hexagons by connecting each point on s1 to its 3 s2 neighbors
    for x, y, s in zip(x_values, y_values, sublattice_values):
        if s == 1:
            for offset in neighbor_offsets:
                neighbor_x = x + offset[0]
                neighbor_y = y + offset[1]

                plt.plot([x, neighbor_x], [y, neighbor_y], "k-", lw=0.5)

            plt.plot(x, y, marker="o", color="k", markersize=1)

        elif s == 2:
            plt.plot(x, y, marker="o", color="w", markersize=1)


def make_contour_plot(
    X, Y, Z, x_values, y_values, sublattice_values, draw_lattice=False
):
    Z_no_nans = np.unique(Z[~np.isnan(Z)])
    Z_no_nans_abs_max = np.abs(Z_no_nans).max()
    Z_no_nans_abs_max_rescaled = Z_no_nans_abs_max / 10

    levels = np.linspace(Z_no_nans.min(), Z_no_nans_abs_max_rescaled, 50)

    plt.contourf(X, Y, Z, levels=levels, cmap="plasma")

    if draw_lattice:
        draw_hexagonal_lattice(x_values, y_values, sublattice_values)

    plt.colorbar()

    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.title("2D Contour Plot")

    plt.show()


def plot_2d_fos(s1=1, lam=0.1, Ef=0.3, draw_lattice=False):
    """
    Use sector method to avoid duplication of equivalent calculations
    and store new calculation results
    """

    ldos_dictionary = {}

    x = []
    y = []
    z = []
    sublattice_tracker = []

    r = 20

    for m in range(-r, r):
        for n in range(-r, r):
            # for s2 in [1, 2]:
            for s2 in [1]:
                m_reduced, n_reduced, s1_reduced, s2_reduced = graphene_sector_method(
                    m, n, s1, s2
                )
                vector_str = "_".join(
                    [str(m_reduced), str(n_reduced), str(s1_reduced), str(s2_reduced)]
                )

                if vector_str not in ldos_dictionary:
                    # calculate
                    ldos = local_density_of_states(
                        lam, Ef, m_reduced, n_reduced, s1_reduced, s2_reduced
                    )
                    # store result
                    ldos_dictionary[vector_str] = ldos
                else:
                    # retrieve result
                    ldos = ldos_dictionary[vector_str]

                x_real, y_real = convert_lattice_vectors_to_real_space(m, n, s2)
                x.append(x_real)
                y.append(y_real)
                z.append(ldos)
                sublattice_tracker.append(s2)

    x_values = np.array(x)
    y_values = np.array(y)
    z_values = np.array(z)
    sublattice_values = np.array(sublattice_tracker)

    # Create a grid for X and Y
    xi = np.linspace(x_values.min(), x_values.max(), 100)
    yi = np.linspace(y_values.min(), y_values.max(), 100)
    X, Y = np.meshgrid(xi, yi)

    # Interpolate Z onto the grid
    Z = griddata((x_values, y_values), z_values, (X, Y), method="cubic")
    Z_symlog = np.sign(Z) * np.log1p(np.abs(Z))
    make_contour_plot(
        X, Y, Z_symlog, x_values, y_values, sublattice_values, draw_lattice
    )


if __name__ == "__main__":
    plot_2d_fos(
        draw_lattice=True,
    )
