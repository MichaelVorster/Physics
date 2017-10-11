# This code calculates the structure function for the velocity and magnetic
# field.  The calculation takes into account that a planar shock can be
# present.

# LIMITATION: code currently only works for Cartesian grids

# Author: Michael Vorster
# Last updated: 11 October 2017

import os
import pyPLUTO as pp
from numpy import (
    arccos,
    cos,
    power,
    random,
    sin,
    square,
    sqrt
)
from matplotlib.pyplot import (
    axis,
    legend,
    loglog,
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
# computational domain.  If not, then shift randomly chosen initial
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
def calculate_difference_squared(D, dimensions, indices, d_index, del_qx):
    points = [[0]*dimensions, [0]*dimensions]
    for i in range(0, dimensions):
        points[0][i] = indices[i]
        points[1][i] = indices[i] + d_index[i]

    B_field = [[0.]*dimensions, [0.]*dimensions]
    points_del_qx = [[0.]*dimensions, [0.]*dimensions]
    for j in range(0, 2):
        for i in range(0, dimensions):
            B_field[j][i] = getattr(D, 'bx' + str(i+1))[
                points[j][0], points[j][1], points[j][2]
            ]
            points_del_qx[j][i] = del_qx[i][
                points[j][0], points[j][1], points[j][2]
            ]

    avrg_B = [0.]*dimensions
    diff_del_qx = [0.]*dimensions
    dot_product = 0

    for i in range(0, dimensions):
        avrg_B[i] = (B_field[0][i] + B_field[1][i])/2.
        diff_del_qx[i] = points_del_qx[1][i] - points_del_qx[0][i]
        dot_product = dot_product + avrg_B[i]*diff_del_qx[i]

    avrg_B_mag = sqrt(sum(square(avrg_B)))
    diff_del_qx_mag = sqrt(sum(square(diff_del_qx)))
    theta = arccos(dot_product/(avrg_B_mag*diff_del_qx_mag))

    diff_qx_para_squared = (diff_del_qx_mag*cos(theta))**2
    diff_qx_perp_squared = (diff_del_qx_mag*sin(theta))**2

    return [diff_qx_para_squared, diff_qx_perp_squared]


# Calculate the correlation function as a function of grid separation
def calculate_structure_function(
    D,
    num_sample_points,
    max_step_size,
    nx,
    del_qx,
    dimensions
):
    structure_function = [[], []]
    for step_size in range(1, max_step_size + 1):
        diff_squared = [0., 0.]
        avrg_diff_squared = [0., 0.]
        for n in range(0, num_sample_points):
            d_index = select_step_direction(dimensions, step_size)
            indices = select_random_point(dimensions, nx, d_index)
            diff_squared_single_point = calculate_difference_squared(
                    D, dimensions, indices, d_index, del_qx
                )
            for i in [0, 1]:
                diff_squared[i] = \
                    diff_squared[i] + diff_squared_single_point[i]
        for i in [0, 1]:
            avrg_diff_squared[i] = diff_squared[i]/num_sample_points
            structure_function[i].append(avrg_diff_squared[i])
    return structure_function


def plot_structure_function(
    structure_function,
    max_step_size,
    file_time,
    region,
    wdir,
    graph_title=''
):
    graph_dir = wdir + 'structure_functions/'
    if not os.path.exists(graph_dir):
        os.makedirs(graph_dir)

    x = range(1, max_step_size+1)
    x_max = max_step_size + 1
    y_max = 1.2*max(structure_function[1])

    index = 2./3.
    norm_x_pos = 20
    normalisation = structure_function[1][norm_x_pos]/power(norm_x_pos, index)
    fit_line = normalisation*power(x, index)

    # loglog(x, structure_function[0])
    loglog(x, structure_function[1], label=r'$SF_{2, \perp}$')
    loglog(x, fit_line, label=r'$x^{2/3}$')
    xlabel(r'Distance  [grid points]')
    ylabel(r'$SF_2$  [scaled units]')
    title(graph_title.title())
    axis([0, x_max, 0, y_max])
    legend()
    savefig(
        graph_dir +
        graph_title.replace(' ', '_') +
        '_structure_function_' +
        region +
        '_t_' +
        str(file_time)
    )
    show()


if __name__ == '__main__':
    # wdir = '/home/mvorster/PLUTO/Shock_turbulence/Results/Run_15_b/'
    # file_time = 0
    # shock_present = 1
    # shock_direction = 'x1'

    wdir = '/home/mvorster/PLUTO/B_Shock_turbulence/bx0.8_cs1.35/'
    file_time = 0
    shock_present = 1
    shock_direction = 'x2'

    region = 'upstream'
    plot_averages = 0
    num_sample_points = 10000

    if not wdir[-1] == '/':
        wdir = wdir + '/'
    D = pp.pload(file_time, w_dir=wdir + 'output/')

    if shock_present:
        if file_time == 0:
            shock_index = 5  # PLUTO simulation has pressure wall at inner boundary
        else:
            shock_index = 260  # locate_shock(D, shock_direction)
    else:
        shock_index = 0
        region = 'no_shock'

    dimensions, nx, max_step_size = grid_info(
        D,
        shock_index,
        shock_direction,
        region
    )

    fluid_quantities = [r'velocity', r'magnetic field']
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
            D,
            num_sample_points,
            max_step_size,
            nx,
            del_qx,
            dimensions
        )

        plot_structure_function(
          structure_function_qx,
          max_step_size,
          file_time,
          region,
          wdir,
          graph_title=fluid_quantity
        )
