def select_starting_position(dimensions, nx, grid_max):
    coordinates = [0]*dimensions
    for dimension in range(0, dimensions):
        coordinates[dimension] = random.random()*grid_max[dimension]

    return coordinates


def calculate_next_position(
    dimensions,
    starting_position,
    B_field,
    step_length,
    grid_max
):
    out_of_domain = 0
    B_dot_product = 0
    for dimension in range(0, dimensions):
        B_dot_product = B_dot_product + square(B_field[dimension])
    B_magnitude = sqrt(B_dot_product)

    next_position = [0]*dimensions
    for dimension in range(0, dimensions):
        step = B_field[dimension]*step_length/B_magnitude
        next_position[dimension] = starting_position[dimension] + step
        if next_position[dimension] > grid_max[dimension]:
        	out_of_domain = 1

    print(next_position)
    	

    return next_position, out_of_domain


def interpolate_magnetic_field(dimensions, grid, B_component_array, position):
    B_field = [0]*dimensions
    for dimension in range(0, dimensions):
        interpolate_b = RegularGridInterpolator(
            points=grid,
            values=B_component_array[dimension]
        )
        B_field[dimension] = interpolate_b(position)

    return B_field


def construct_field_line(D):
    dimensions = 2
    nx = [len(D.x2), len(D.x3)]
    step_length = (max(D.x1) - min(D.x1))/(len(D.x1) - 1.)/2.
    number_of_steps = 1

    grid = [D.x2, D.x3]
    grid_max = [max(D.x2), max(D.x3)]
    B_component_array = [D.bx2[512, :, :], D.bx3[512, :, :]]

    starting_position = select_starting_position(dimensions, nx, grid_max)
    B_field = interpolate_magnetic_field(
        dimensions,
        grid,
        B_component_array,
        starting_position
    )

    print(starting_position)

    quiver(
        D.x2[60:70],
        D.x3[60:70],
        D.bx2[512, 60:70, 60:70],
        D.bx3[512, 60:70, 60:70],
        units='width'
    )

    out_of_domain = 0
    step = 1
    while not out_of_domain and step <= number_of_steps:
        next_position, out_of_domain = calculate_next_position(
            dimensions,
            starting_position,
            B_field,
            step_length,
            grid_max
        )
        step += 1
        '''
        B_field = interpolate_magnetic_field(
            dimensions,
            grid,
            B_component_array,
            next_position
        )
        starting_position = next_position

        quiver(
            next_position[0],
            next_position[1],
            B_field[0],
            B_field[1],
            color='r'
        )

        show()
        '''

if __name__ == '__main__':
    from matplotlib.pyplot import quiver, show
    from numpy import (
        random,
        square,
        sqrt
    )
    from pyPLUTO import pload
    from scipy.interpolate import RegularGridInterpolator

    wdir = '/home/mvorster/PLUTO/Shock_turbulence/output/'
    D = pload(0, w_dir=wdir)

    construct_field_line(D)
