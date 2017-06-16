# This code calculates the structure function spectrum for the density,
# velocity, and magnetic field.  It is assumed that there is a planar shock
# present, and that this shock propagates in the x-direction

# LIMITATION: code currently only works for Cartesian grids

# Author: Michael Vorster
# Last updated: 16 June 2017

from fitcurve import savitzky_golay
from numpy import (
    argmax,
    array,
    average,
    random
)
from matplotlib.pyplot import (
    axis,
    plot,
    savefig,
    show,
    title,
    xlabel,
    ylabel
)


def grid_info(D):
    # Number of dimensions
    dimensions = 3
    nx = (len(D.x1), len(D.x2), len(D.x3))
    dx = [
        (max(D.x1) - min(D.x1))/(nx[0] - 1.),
        (max(D.x2) - min(D.x2))/(nx[1] - 1.),
        (max(D.x3) - min(D.x3))/(nx[2] - 1.)
    ]
    # we will only sample points that are not more than half a grid size
    # separated.  If grid does not have an equal size in all directions,
    # then the step size is determined by the smallest grid direction
    max_step_size = min(nx)/2

    return dimensions, nx, dx, max_step_size


# Select step direction
def select_step_direction(dimensions, step_size):
    # step should be only either in x1, x2, or x3 direction
    # step can be either in positive or negative direction
    index = 0
    while index == 0:
        index = random.random_integers(-1*dimensions, 1*dimensions)
    d_index = [0]*dimensions
    if index < 0:
        d_index[abs(index)-1] = -1*step_size
    else:
        d_index[index-1] = 1*step_size

    return d_index


# Calculate the average profile of a fluid quantity along the x-axis, i.e.,
# the fluid quantity is averaged along the y- and z-axes.  In general the
# averaged profile will not be smooth.  Therefore, it is necessary to fit a
# smooth function to the profile.  However, the presence of a shock complicates
# the problem.
#
# The approach is to fit the shocked region with a smooth curve.  Upstream of
# the shock the fluid quantity is again averaged, but this time over all axes.
# The second averaging should ensure a smooth profile upstream of the shock.
def average_fluid_field(x1, qx, nx, fluid_quantity):
    avrg_qx = [0]*nx[0]
    for i in range(0, nx[0]):
        avrg_qx[i] = average(qx[i, :, :])

    original_profile = list(avrg_qx)

    # it is usefule to fit the data a few grid points beyond shock position
    shock_position = argmax(avrg_qx) + 1
    downstream_fit = savitzky_golay(
        avrg_qx[0:shock_position], 41, 1
    )
    upstream_fit = average(qx[shock_position+1:, :, :])

    for i in range(0, shock_position):
        avrg_qx[i] = downstream_fit[i]
    for i in range(shock_position, nx[0]):
        avrg_qx[i] = upstream_fit

    y_max = 1.2*max(avrg_qx)
    plot(x1, avrg_qx, 'k', x1, original_profile, 'r')
    xlabel(r'x')
    ylabel(r'Fitted (black), Original (red)')
    title('Average ' + fluid_quantity)
    axis([0, 1, -0.5, y_max])
    savefig(fluid_quantity+'_x_average')
    show()

    return avrg_qx


# Select random point
def select_random_point(dimensions, nx, d_index):
    indices = [0]*dimensions
    for dimension in range(0, dimensions):
        indices[dimension] = random.random_integers(0, nx[dimension]-1)

    test_index_step(dimensions, nx, d_index, indices)

    return indices


# Test that second point chosen based on step size is still inside
# computational domain.  If not, then move shift randomly chosen initial
# coordinate.
def test_index_step(dimensions, nx, d_index, indices):
    for dimension in range(0, dimensions):
        if indices[dimension] + d_index[dimension] > (nx[dimension] - 1):
            new_index = random.random_integers(
                0, nx[dimension] - d_index[dimension] - 1
            )
            indices[dimension] = new_index
    return indices


# Calculate the square of the difference of a fluid quantity between
# two selected points.  If fluid quantity is a vector, then the difference
# is calculated for every vector component
def calculate_difference_squared(dimensions, indices, d_index, del_qx):
    difference = [0]*len(del_qx)
    difference_squared = 0
    for i in range(0, len(del_qx)):
        x11 = indices[0]
        x12 = indices[0] + d_index[0]
        x21 = indices[1]
        x22 = indices[1] + d_index[1]
        x31 = indices[2]
        x32 = indices[2] + d_index[2]

        difference[i] = del_qx[i][x11, x21, x31] - del_qx[i][x12, x22, x32]
        difference_squared = difference_squared + difference[i]**2
    return difference_squared


# Calculate the correlation function as a function of grid separation
def calculate_structure_function(
    num_sample_points,
    nx,
    dimensions,
    max_step_size,
    del_qx
):
    structure_function = []
    for step_size in range(1, max_step_size + 1):
        diff_squared = 0
        for n in range(1, num_sample_points):
            d_index = select_step_direction(dimensions, step_size)
            indices = select_random_point(dimensions, nx, d_index)
            diff_squared_single_point = calculate_difference_squared(
                dimensions, indices, d_index, del_qx
            )
            diff_squared = diff_squared + diff_squared_single_point
            avrg_diff_squared = diff_squared/num_sample_points
        structure_function.append(avrg_diff_squared)
    return structure_function


def plot_structure_function(
    structure_function,
    max_step_size,
    graph_title=''
):
    x = range(1, max_step_size+1)
    x_max = max_step_size + 1
    y_max = 1.2*max(structure_function)

    plot(x, structure_function)
    xlabel(r'x')
    ylabel(r'$S(x)$')
    title(graph_title + ' structure function')
    axis([0, x_max, 0, y_max])
    savefig(graph_title+'_structure_function')
    show()


def density_turbulence_spectrum(D, num_sample_points):
    dimensions, nx, dx, max_step_size = grid_info(D)

    # Average values
    avrg_rho = average_fluid_field(
        D.x1,
        D.rho,
        nx,
        'Density'
    )

    # Turbulent variations
    del_rho = array(list(D.rho))
    for i in range(0, nx[0]):
        del_rho[i, :, :] = array(D.rho[i, :, :] - avrg_rho[i])

    del_rho = [del_rho]

    structure_function_rho = calculate_structure_function(
        num_sample_points,
        nx,
        dimensions,
        max_step_size,
        del_rho
    )

    plot_structure_function(
      structure_function_rho,
      max_step_size,
      graph_title=r'Density'
    )


def velocity_turbulence_spectrum(D, num_sample_points):
    dimensions, nx, dx, max_step_size = grid_info(D)

    # Average values
    avrg_vx1 = average_fluid_field(
        D.x1,
        D.vx1,
        nx,
        '$V_{x_1}$'
    )
    avrg_vx2 = average(D.vx2)
    avrg_vx3 = average(D.vx3)

    # Turbulent variations
    del_vx1 = array(list(D.vx1))
    for i in range(0, nx[0]):
        del_vx1[i, :, :] = array(D.vx1[i, :, :] - avrg_vx1[i])
    del_vx2 = array(D.vx2)-avrg_vx2
    del_vx3 = array(D.vx3)-avrg_vx3

    del_vx = [
        del_vx1,
        del_vx2,
        del_vx3
    ]

    structure_function_vx = calculate_structure_function(
        num_sample_points,
        nx,
        dimensions,
        max_step_size,
        del_vx
    )

    plot_structure_function(
      structure_function_vx,
      max_step_size,
      graph_title=r'Velocity'
    )


if __name__ == "__main__":
    import os
    import pyPLUTO as pp

    plutodir = os.environ['PLUTO_DIR']
    wdir = '/home/mvorster/PLUTO/Shock_turbulence/Results/Run_8/output/'

    D = pp.pload(7, w_dir=wdir)
    num_sample_points = 10000

    density_turbulence_spectrum(D, num_sample_points)
    velocity_turbulence_spectrum(D, num_sample_points)
