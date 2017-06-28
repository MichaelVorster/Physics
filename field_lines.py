def select_starting_indices(dimensions, nx):
    indices = [0]*dimensions
    for dimension in range(0, dimensions):
        indices[dimension] = random.random_integers(0, nx[dimension]-1)

    return indices


def calculate_next_position(D, indices, step_length):
    ni, nj, nk = indices

    B_magnitude = sqrt(
        square(D.bx1[ni, nj, nk]) +
        square(D.bx2[ni, nj, nk]) +
        square(D.bx3[ni, nj, nk])
    )
    
    step = array([
        D.bx1[ni, nj, nk],
        D.bx2[ni, nj, nk],
        D.bx3[ni, nj, nk]
    ])*step_length/B_magnitude
    
    next_position = [
        D.x1[ni] + step[0],
        D.x2[nj] + step[1],
        D.x3[nk] + step[2]
    ]    

    return next_position


# The next point should be somewhere in a cube defined by eight grid points
def calculate_surrounding_indices(next_position, step_length):
    ni = floor(next_position[0]/step_length)
    nj = floor(next_position[1]/step_length)
    nk = floor(next_position[2]/step_length)


def factory_interpolate(D, bx):
    return RegularGridInterpolator(
        points=[
            D.x1, D.x2, D.x3
        ],
        values=bx
    )


def interpolate_magnetic_field(D, position, step_length):
    interpolate_bx1 = factory_interpolate(D, D.bx1)
    interpolate_bx2 = factory_interpolate(D, D.bx2)
    interpolate_bx3 = factory_interpolate(D, D.bx3)
 
    bx1 = interpolate_bx1(position)
    bx2 = interpolate_bx2(position)
    bx3 = interpolate_bx3(position)

    return [bx1, bx2, bx3]


def construct_field_line(D, dimensions, nx):
    step_length = (max(D.x1) - min(D.x1))/(len(D.x1) - 1.)/2.

    starting_indices = select_starting_indices(dimensions, nx)
    next_position = calculate_next_position(D, starting_indices, step_length)
    next_indices = calculate_surrounding_indices(next_position, step_length)

    B_field = interpolate_magnetic_field(D, next_position, step_length)
    
    nx, ny, nz = starting_indices

    print('Starting point: \n')
    print(D.x1[nx], D.x2[ny], D.x3[nz])
    print('\n')

    print('Starting point B field: \n')
    print(D.bx1[nx, ny, nz], D.bx2[nx, ny, nz], D.bx3[nx, ny, nz])
    print('\n')

    print('New point: \n')
    print(next_position)
    print('\n')

    print('New point B field: \n')
    print(B_field)
    print('\n')


if __name__ == '__main__':

    from math import floor
    from matplotlib.pyplot import *
    from numpy import (
        array,
        random,
        square,
        sqrt
    )
    from pyPLUTO import pload
    from scipy.interpolate import RegularGridInterpolator


    wdir = '/home/mvorster/PLUTO/Shock_turbulence/output/'
    D = pload(0 ,w_dir=wdir)

    dimensions = 3
    nx = [len(D.x1), len(D.x2), len(D.x3)]

    construct_field_line(D, dimensions, nx)