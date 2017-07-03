# Code to investigate Richardson diffusion of magnetic field lines.
# Two field lines are constructed starting from two points separated by
# some initial distance.  The separation between the end points of the
# field lines are then calculated.  This procedure is repeated for a large
# number of field line pairs
#
# AUTHOR: Michael Vorster
#
# LAST UPDATED: 03 July 2017


def select_starting_positions(dimensions, nx, grid_min, grid_max, separation):
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
        point_1_in_domain = False
        point_2_in_domain = False
        while not (point_1_in_domain and point_2_in_domain):
            starting_positions[0][dimension] = \
                random.random()*grid_max[dimension]
            starting_positions[1][dimension] = \
                starting_positions[0][dimension] + \
                dx[dimension]

            point_1_in_domain, point_2_in_domain = point_in_domain(
                [
                    starting_positions[0][dimension],
                    starting_positions[1][dimension]
                ],
                grid_min[dimension],
                grid_max[dimension]
            )

    return starting_positions


def calculate_next_position(
    dimensions,
    starting_position,
    B_field,
    step_length,
    grid_min,
    grid_max
):
    out_of_domain = 0
    B_dot_product = 0
    for dimension in range(0, dimensions):
        B_dot_product = B_dot_product + square(B_field[dimension][0])
    B_magnitude = sqrt(B_dot_product)

    next_position = [0.]*dimensions
    for dimension in range(0, dimensions):
        # B_field[dimension][0] - B components are arrays inside an array
        step = B_field[dimension][0]*step_length/B_magnitude
        next_position[dimension] = starting_position[dimension] + step

        if next_position[dimension] < grid_min[dimension] or \
            next_position[dimension] > grid_max[dimension]:
                out_of_domain = 1

    return next_position, out_of_domain


def point_in_domain(positions, grid_min, grid_max):
    in_domain = []
    for position in positions:
        in_domain.append(position > grid_min and position < grid_max)

    return in_domain


def interpolate_magnetic_field(dimensions, grid, B_component_array, position):
    B_field = [0.]*dimensions
    for dimension in range(0, dimensions):
        interpolate_b = RegularGridInterpolator(
            points=grid,
            values=B_component_array[dimension]
        )
        B_field[dimension] = interpolate_b(position)

    return B_field


def remove_additional_point(
    dimensions,
    field_line_coordinates,
    field_line_components
):
    # if second point is out of domain, remove the last entry from
    # first point.  Otherwise point one will have one more entry
    for dimension in range(0, dimensions):
        field_line_coordinates[0][dimension].pop()
        field_line_components[0][dimension].pop()


def get_grid_info(dimensions, D):
    # assumes that the grid separation is the same in all directions
    step_length = (max(D.x1) - min(D.x1))/(len(D.x1) - 1.)/2.

    # for plotting purposes arrays must have same dimensions.
    # Only x2 and x3 have the same dimensions
    if dimensions == 2:
        x_slice = len(D.x1)/2
        grid = [D.x2, D.x3]
        B_component_array = [D.bx2[x_slice, :, :], D.bx3[x_slice, :, :]]
    if dimensions == 3:
        grid = [D.x1, D.x2, D.x3]
        B_component_array = [D.bx1, D.bx2, D.bx3]

    nx = [0]*dimensions
    grid_min = [0.]*dimensions
    grid_max = [0.]*dimensions
    for dimension in range(0, dimensions):
        nx[dimension] = len(grid[dimension])
        grid_min[dimension] = min(grid[dimension])
        grid_max[dimension] = max(grid[dimension])

    return nx, step_length, grid, grid_min, grid_max, B_component_array


def construct_field_lines(dimensions, D, separation, number_of_steps):
    nx, step_length, grid, grid_min, grid_max, B_component_array = \
        get_grid_info(dimensions, D)

    starting_positions = select_starting_positions(
        dimensions,
        nx,
        grid_min,
        grid_max,
        separation
    )

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
    while step <= number_of_steps and not out_of_domain:
        for point in (0, 1):
            next_position, out_of_domain = calculate_next_position(
                dimensions,
                starting_positions[point],
                B_field[point],
                step_length,
                grid_min,
                grid_max
            )

            if out_of_domain:
                if step == 1:
                    # pair of points should be discarded if points are
                    # already out of domain on first step
                    step = -1
                if point == 1:
                    remove_additional_point(
                        dimensions,
                        field_line_coordinates,
                        field_line_components
                    )
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

    return field_line_coordinates, field_line_components, step


def calculate_field_line_separation(dimensions, field_lines):
    diffusion = [[], []]
    number_of_steps = len(field_lines[0][0])
    for step in range(0, number_of_steps):
        separation_vector = array([0.]*dimensions)
        for dimension in range(0, dimensions):
            separation_vector[dimension] = \
                field_lines[1][dimension][step] - \
                field_lines[0][dimension][step]
        separation_magnitude = sqrt(sum(square(separation_vector)))
        diffusion[0].append(step)
        diffusion[1].append(separation_magnitude)

    return diffusion


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

    show()


def add_result(number_of_steps, average_diffusion, diffusion):
    number_new_steps = len(diffusion[0])
    for step in range(0, number_new_steps):
        average_diffusion[1][step] = \
            average_diffusion[1][step] + diffusion[1][step]

    return average_diffusion


def plot_diffusion(number_of_separations, average_diffusion_per_separation):
    Richardson_x = range(0, 100)
    Richardson_y = power(Richardson_x, 0.5)
    for separation in range(0, number_of_separations):
        loglog(
            average_diffusion_per_separation[separation][0],
            average_diffusion_per_separation[separation][1]
        )
    loglog(Richardson_x, 1e-3*Richardson_y)
    show()


def calculate_B_field_line_diffusion(D):
    dimensions = 3
    number_of_pairs = 30
    number_of_steps = 20  # along B
    number_of_separations = 2
    inital_separations = linspace(
        0.001,
        0.004,
        number_of_separations
    )

    average_diffusion_per_separation = []
    for separation in inital_separations:
        number_usable_pairs = 0
        average_diffusion = [
            range(0, number_of_steps),
            [0]*number_of_steps
        ]
        for pair in range(0, number_of_pairs):
            field_line_coordinates, field_line_components, step = \
                construct_field_lines(
                    dimensions,
                    D,
                    separation,
                    number_of_steps
                )
            if step > -1:
                diffusion = calculate_field_line_separation(
                    dimensions,
                    field_line_coordinates
                )
                add_result(number_of_steps, average_diffusion, diffusion)
                number_usable_pairs += 1

            # plot field lines
            max_plot_number = -1
            plot_number = 0
            if dimensions == 2 and plot_number <= max_plot_number:
                plot_field_lines(
                    D,
                    field_line_coordinates,
                    field_line_components
                )
            plot_number += 1

        average_diffusion_per_separation.append([
            average_diffusion[0],
            array(average_diffusion[1])/number_usable_pairs
        ])

    plot_diffusion(number_of_separations, average_diffusion_per_separation)


if __name__ == '__main__':
    from math import (
        sin,
        cos
    )
    from matplotlib.pyplot import (
        axis,
        loglog,
        quiver,
        show
    )
    from numpy import (
        array,
        linspace,
        power,
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
