# Code to investigate Richardson diffusion of magnetic field lines.
# Two field lines are constructed starting from two points separated by
# some initial distance.  The separation between the end points of the
# field lines are then calculated.  This procedure is repeated for a large
# number of field line pairs
#
# AUTHOR: Michael Vorster
#
# LAST UPDATED: 30 June 2017


def select_starting_positions(dimensions, nx, grid_max, separation):
    starting_positions = [
        [0.]*dimensions,
        [0.]*dimensions
    ]

    phi = random.random()*3.14159
    if dimensions == 2:
        theta = 3.14159/2.
    if dimensions == 3:
        theta = random.random()*3.14159

    dx = [
        separation*cos(phi)*sin(theta),
        separation*sin(phi)*sin(theta),
        separation*cos(theta)
    ]

    for dimension in range(0, dimensions):
        in_domain = 0
        while not in_domain:
            starting_positions[0][dimension] = \
                random.random()*grid_max[dimension]
            starting_positions[1][dimension] = \
                dx[dimension] + \
                starting_positions[0][dimension]
            if starting_positions[1][dimension] < grid_max[dimension]:
                in_domain = 1

    return starting_positions


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
        B_dot_product = B_dot_product + square(B_field[dimension][0])
    B_magnitude = sqrt(B_dot_product)

    next_position = [0.]*dimensions
    for dimension in range(0, dimensions):
        # 'B_field[dimension][0]' - this is a data type conversion of sorts
        # B_field components are arrays inside an array
        step = B_field[dimension][0]*step_length/B_magnitude
        next_position[dimension] = starting_position[dimension] + step
        if next_position[dimension] > grid_max[dimension]:
            out_of_domain = 1

    return next_position, out_of_domain


def interpolate_magnetic_field(dimensions, grid, B_component_array, position):
    B_field = [0.]*dimensions
    for dimension in range(0, dimensions):
        interpolate_b = RegularGridInterpolator(
            points=grid,
            values=B_component_array[dimension]
        )
        B_field[dimension] = interpolate_b(position)

    return B_field


def construct_field_lines(dimensions, D, separation):
    number_of_steps = 1

    # assumes that the grid separation is the same in all directions
    step_length = (max(D.x1) - min(D.x1))/(len(D.x1) - 1.)/2.
    x_slice = len(D.x1)/2

    # for plotting purposes arrays must have same dimensions.
    # Only x2 and x3 have the same dimensions
    if dimensions == 2:
        grid = [D.x2, D.x3]
        B_component_array = [D.bx2[x_slice, :, :], D.bx3[x_slice, :, :]]
    if dimensions == 3:
        grid = [D.x1, D.x2, D.x3]
        B_component_array = [D.bx1, D.bx2, D.bx3]

    nx = [0]*dimensions
    grid_max = [0.]*dimensions
    for dimension in range(0, dimensions):
        nx[dimension] = len(grid[dimension])
        grid_max[dimension] = max(grid[dimension])

    starting_positions = select_starting_positions(
        dimensions,
        nx,
        grid_max,
        separation
    )
    '''
    x1 = 0.060
    y1 = 0.060
    x2 = x1 + separation
    y2 = y1
    starting_positions = [[x1, y1], [x2, y2]]
    '''
    # starting_positions = [[0.068, 0.05908203125], [0.068, 0.05908203125]]

    B_field = [[], []]
    for i in (0, 1):
        B_field[i] = interpolate_magnetic_field(
            dimensions,
            grid,
            B_component_array,
            starting_positions[i]
        )

    field_line_coordinates = [[], []]
    field_line_components = [[], []]
    for point in (0, 1):
        for dimension in range(0, dimensions):
            field_line_coordinates[point].append([])
            field_line_components[point].append([])

    out_of_domain = 0
    step = 1
    while step <= number_of_steps:
        for point in (0, 1):
            next_position, out_of_domain = calculate_next_position(
                dimensions,
                starting_positions[point],
                B_field[point],
                step_length,
                grid_max
            )
            if out_of_domain:
                break 
            B_field[point] = interpolate_magnetic_field(
                dimensions,
                grid,
                B_component_array,
                next_position
            )
            starting_positions[point] = next_position[:]

            for dimension in range(0, dimensions):
                field_line_coordinates[point][dimension].append(
                    starting_positions[point][dimension]
                )
                field_line_components[point][dimension].append(
                    B_field[point][dimension]
                )

        step += 1

    return field_line_coordinates, field_line_components


def calculate_field_line_separation(dimensions, field_lines):
    separation_vector = array([0.]*dimensions)
    for dimension in range(0, dimensions):
        separation_vector[dimension] = \
            field_lines[1][dimension][-1] - field_lines[0][dimension][-1]

    separation_magnitude = sqrt(sum(square(separation_vector)))

    return separation_magnitude


def plot_field_lines(
    D,
    field_line_coordinates,
    field_line_components
):
    begin = 0
    end = 255
    x_slice = 512

    axis('equal')
    quiver(
        D.x2[begin:end],
        D.x3[begin:end],
        D.bx2[x_slice, begin:end, begin:end].T,
        D.bx3[x_slice, begin:end, begin:end].T,
        units='xy'
    )

    colour = ['g', 'b']
    for i in (0, 1):
        axis('equal')
        quiver(
            field_line_coordinates[i][0],
            field_line_coordinates[i][1],
            field_line_components[i][0],
            field_line_components[i][1],
            units='xy',
            color=colour[i]
        )

    # show()


def calculate_B_field_line_diffusion(D):
    dimensions = 2
    number_of_pairs = 3
    inital_separations = linspace(
        0.004,
        0.02,
        3
    )
    final_separations = []

    for separation in inital_separations:
        for pair in range(0, number_of_pairs):
            field_line_coordinates, field_line_components = construct_field_lines(  # noqa
                dimensions,
                D,
                separation
            )
            separation = calculate_field_line_separation(
                dimensions,
                field_line_coordinates
            )
            final_separations.append(separation)

            # print(separation)

            # starting_positions = select_starting_positions(dimensions, nx, grid_max, separation)  # noqa
            # starting_positions = [[0.068, 0.05908203125], [0.060, 0.05908203125]]  # noqa
            max_plot_number = 2
            plot_number = 0
            if dimensions == 2 and plot_number <= max_plot_number:
                plot_field_lines(
                    D,
                    field_line_coordinates,
                    field_line_components
                )
            plot_number += 1


if __name__ == '__main__':
    from math import (
        sin,
        cos
    )
    from matplotlib.pyplot import (
        axis,
        quiver,
        show
    )
    from numpy import (
        array,
        linspace,
        random,
        square,
        sqrt,
        sum
    )
    from pyPLUTO import pload
    from scipy.interpolate import RegularGridInterpolator

    wdir = '/home/mvorster/PLUTO/Shock_turbulence/output/'
    D = pload(0, w_dir=wdir)

    calculate_B_field_line_diffusion(D)
