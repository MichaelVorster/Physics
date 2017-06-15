# This code calculates the structure function spectrum for the density,
# velocity, and magnetic field.

# LIMITATION: code currently only works for Cartesian grids

# Author: Michael Vorster
# Last updated: 13 June 2017

from numpy import *
from matplotlib.pyplot import *


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


# Calculate correlation between selected points
def calculate_correlation_single_point(dimensions, indices, d_index, del_qx):
    vector_diff = [0]*len(del_qx)
    correlation = 0
    for i in range(0, len(del_qx)):

        x1 = indices[0]
        x2 = indices[0] + d_index[0]
        y1 = indices[1]
        y2 = indices[1] + d_index[1]
        z1 = indices[2]
        z2 = indices[2] + d_index[2]

        vector_diff[i] = del_qx[i][x1,y1,z1] - del_qx[i][x2,y2,z2]
        correlation = correlation + vector_diff[i]**2
    return correlation


# Calculate the correlation function as a function of grid separation
def calculate_correlation_vector(
    num_sample_points,
    nx,
    dimensions,
    max_step_size,
    del_qx
):
    correlation_vector = []
    for step_size in range(1, max_step_size + 1):
        correlation = 0
        for n in range(1, num_sample_points):
            d_index = select_step_direction(dimensions, step_size)
            indices = select_random_point(dimensions, nx, d_index)
            correlation_single_point = calculate_correlation_single_point(
                dimensions, indices, d_index, del_qx
            )
            correlation = correlation + correlation_single_point
        avrg_correlation = correlation/num_sample_points
        correlation_vector.append(avrg_correlation)
    return correlation_vector


def plot_function(correlation_vector, max_step_size, graph_title=''):
    x = range(1, max_step_size+1)
    x_max = max_step_size + 1
    y_max = 1.2*max(correlation_vector)

    plot(x, correlation_vector)
    xlabel(r'x')
    ylabel(r'$S(x)$')
    title(graph_title)
    axis([0,x_max,0,y_max])
    show() 


def density_turbulence_spectrum(D, num_sample_points):
    dimensions, nx, dx, max_step_size = grid_info(D)
    # Average values
    avrg_rho = average(D.rho)
    # Turbulent variations
    del_rho = [
        array(D.rho)-avrg_rho 
    ]

    correlation_vector_rho = calculate_correlation_vector(
        num_sample_points,
        nx,
        dimensions,
        max_step_size,
        del_rho
    )

    plot_function(
      correlation_vector_rho,
      max_step_size,
      graph_title=r'Density'
    )


def velocity_turbulence_spectrum(D, num_sample_points):
    dimensions, nx, dx, max_step_size = grid_info(D)
    # Average values
    avrg_vx1 = average(D.vx1)
    avrg_vx2 = average(D.vx2)
    avrg_vx3 = average(D.vx3)
    # Turbulent variations
    del_vx = [
        array(D.vx1)-avrg_vx1,
        array(D.vx2)-avrg_vx2,
        array(D.vx3)-avrg_vx3,
    ]
    correlation_vector_vx = calculate_correlation_vector(
        num_sample_points,
        nx,
        dimensions,
        max_step_size,
        del_vx
    )

    plot_function(
      correlation_vector_vx,
      max_step_size,
      graph_title=r'Velocity'
    )


if __name__ == "__main__":
    import os
    import pyPLUTO as pp

    plutodir = os.environ['PLUTO_DIR']
    wdir = '/home/mvorster/PLUTO/Shock_turbulence/output/'

    D = pp.pload(5, w_dir=wdir)
    num_sample_points = 10000

    density_turbulence_spectrum(D, num_sample_points)
    velocity_turbulence_spectrum(D, num_sample_points)
