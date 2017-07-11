# This code calculates the structure function spectrum for the density,
# velocity, and magnetic field.  It is assumed that there is a planar shock
# present, and that this shock propagates in the x-direction

# LIMITATION: code currently only works for Cartesian grids

# Author: Michael Vorster
# Last updated: 16 June 2017

from numpy import (
    argmax,
    array,
    average,
    power,
    random
)
from matplotlib.pyplot import (
    axis,
    loglog,
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
    flag_avrg,
    graph_title=''
):
    x = range(1, max_step_size+1)
    len_x = len(x)
    x_max = max_step_size + 1
    y_max = 1.2*max(structure_function)

    norm = 1./(x[len_x-1]**(2./3.))*structure_function[len_x-1]*2
    power_law = power(x, 2./3.)*norm

    loglog(x, structure_function, 'k', x, power_law, 'b')
    xlabel(r'x')
    ylabel(r'$S(x)$')
    title(graph_title + ' structure function')
    axis([0, x_max, 0, y_max])
    savefig(
        'Initial/Initial_'+graph_title+'_structure_function_'+str(flag_avrg)
    )
    show()


def density_turbulence_spectrum(D, num_sample_points, flag_avrg):
    dimensions, nx, dx, max_step_size = grid_info(D)

    # Turbulent variations
    avrg_rho = average(D.rho)*flag_avrg
    del_rho = [array(D.rho) - avrg_rho]

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
      flag_avrg,
      graph_title=r'Density'
    )


def velocity_turbulence_spectrum(D, num_sample_points, flag_avrg):
    dimensions, nx, dx, max_step_size = grid_info(D)

    # Average values
    avrg_vx1 = average(D.vx1)*flag_avrg
    avrg_vx2 = average(D.vx2)*flag_avrg
    avrg_vx3 = average(D.vx3)*flag_avrg

    # Turbulent variations
    del_vx1 = array(D.vx1) - avrg_vx1
    del_vx2 = array(D.vx2) - avrg_vx2
    del_vx3 = array(D.vx3) - avrg_vx3

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
      flag_avrg,
      graph_title=r'Velocity'
    )


if __name__ == "__main__":
    import os
    import pyPLUTO as pp

    plutodir = os.environ['PLUTO_DIR']
    wdir = '/home/mvorster/PLUTO/Shock_turbulence/Results/Run_8/output/'

    D = pp.pload(0, w_dir=wdir)
    num_sample_points = 10000

    for flag_avrg in (0, 1): 
        density_turbulence_spectrum(D, num_sample_points, flag_avrg)
        velocity_turbulence_spectrum(D, num_sample_points, flag_avrg)
