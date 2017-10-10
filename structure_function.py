# This code calculates the structure function spectrum for the density,
# velocity, and magnetic field.  It is assumed that there is a planar shock
# present, and that this shock propagates in the x-direction

# LIMITATION: code currently only works for Cartesian grids

# Author: Michael Vorster
# Last updated: 16 June 2017

import os
import pyPLUTO as pp
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
from shock_tools import (
    locate_shock,
    remove_average_fluid_component
)


def grid_info(
    D,
    shock_index,
    shock_direction,
    region
):
    offset = 10
    dimensions = 3

    if region == 'upstream':
        if shock_direction == 'x1':
            grid = [D.x1[shock_index + offset:], D.x2, D.x3]
        elif shock_direction == 'x2':
            grid = [D.x1, D.x2[shock_index + offset:], D.x3]
        else:
            grid = [D.x1, D.x2, D.x3[shock_index + offset:]]    
    elif region == 'downstream':
        if shock_direction == 'x1':
            grid = [D.x1[:shock_index - offset + 1], D.x2, D.x3]
        elif shock_direction == 'x2':
            grid = [D.x1, D.x2[:shock_index - offset + 1], D.x3]
        else:
            grid = [D.x1, D.x2, D.x3[:shock_index - offset + 1]]
    else:
        grid = [D.x1, D.x2, D.x3]

    nx = (len(grid[0]), len(grid[1]), len(grid[2]))
    # we will only sample points that are separated by a distance that is
    # equal to (or less) than half the size of the region under consideration
    max_step_size = min(nx)/2

    return dimensions, nx, max_step_size


def limit_quantity_to_region(
    del_qx,
    region,
    nx
):
    for i in range(len(del_qx)):
        if region == 'upstream':
            del_qx[i] = del_qx[i][-nx[0]:, -nx[1]:, -nx[2]:]
        elif region == 'downstream':
            del_qx[i] = del_qx[i][0:nx[0], 0:nx[1], 0:nx[2]]

    return del_qx


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
    max_step_size,
    nx,
    del_qx,
    dimensions
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


if __name__ == '__main__':
    wdir = '/home/mvorster/PLUTO/Shock_turbulence/Results/Run_8/output/'
    file_time = 7
    shock_present = 1
    shock_direction = 'x1'
    region = 'upstream'
    plot_averages = 0

    if shock_present:
        shock_index = locate_shock(D)
    else:
        shock_index = 0
        region = 'no_shock'

    D = pp.pload(file_time, w_dir=wdir)
    num_sample_points = 10000

    dimensions, nx, max_step_size = grid_info(
        D,
        shock_index,
        shock_direction,
        region
    )

    fluid_quantities = [r'density', r'velocity', r'magnetic field']
    fluid_quantities = [r'velocity']
    for fluid_quantity in fluid_quantities:
        del_qx = remove_average_fluid_component(
            D,
            fluid_quantity,
            shock_index,
            shock_direction,
            wdir,
            plot_averages
        )

        if shock_present:
            del_qx = limit_quantity_to_region(
                del_qx,
                region,
                nx
            )

        structure_function_qx = calculate_structure_function(
            num_sample_points,
            max_step_size,
            nx,
            del_qx,
            dimensions
        )

        plot_structure_function(
          structure_function_qx,
          max_step_size,
          graph_title=fluid_quantity.title()
        )
