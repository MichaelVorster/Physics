# Code to investigate diffusion of magnetic field lines.  Two field lines are
# constructed starting from two points separated by some initial distance.
# The separation between the field lines are calculated at set intervals along
# the lines.  This procedure is repeated for a number of field line pairs.
#
# AUTHOR: Michael Vorster
#
# LAST UPDATED: 09 October 2017


# If shock is along x-axis, the following changes are needed 
# (might be incomplete):
# Line 127 and 132
# In 'get_grid_info', especially the array slices


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
    arccos,
    array,
    mod,
    random,
    square,
    sqrt
)
from pyPLUTO import pload
from scipy.interpolate import RegularGridInterpolator
from shock_tools import locate_shock


def select_starting_positions(dimensions, nx, grid_min, grid_max, separation):
    starting_positions = [
        [0.]*dimensions,
        [0.]*dimensions
    ]

    phi = random.random()*3.14159
    if dimensions == 2:
        theta = 3.14159/2.
    if dimensions == 3:
        theta = random.random()*3.14159/2.

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
    grid_max,
    shock_index
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

    # periodic boundaries in x and z, and in y when no shock is present
    if shock_index > 0:
        if next_position[1] < grid_min[1] or next_position[1] > grid_max[1]:
            out_of_domain = 1

    return next_position, out_of_domain


def point_in_domain(positions, grid_min, grid_max):
    in_domain = []
    for position in positions:
        in_domain.append(position > grid_min and position < grid_max)

    return in_domain


def interpolate_magnetic_field(
    dimensions,
    grid,
    B_component_array,
    position,
    shock_index
):
    B_field = [0.]*dimensions
    position_temp = [0.]*dimensions

    # periodic boundaries along x and z (and y if no shock is present)
    dimensions_range = [0, 2]
    if shock_index == 0 or dimensions == 2:
        dimensions_range = range(0, dimensions)

    # periodic boundaries along y and z (and x if no shock is present)
    position_temp[1] = position[1]
    for dimension in dimensions_range:
        max_grid = max(grid[dimension])
        min_grid = min(grid[dimension])
        grid_extension = max_grid - min_grid
        if position[dimension] > max_grid:
            correction = mod(
                abs(position[dimension] - max_grid), grid_extension
            )
            position_temp[dimension] = min_grid + correction
        elif position[dimension] < min_grid:
            # inner boundary is not at zero
            correction = mod(
                abs(position[dimension] - min_grid), grid_extension
            )
            position_temp[dimension] = max_grid - correction
        else:
            position_temp[dimension] = position[dimension]

    for dimension in range(0, dimensions):
        interpolate_b = RegularGridInterpolator(
            points=grid,
            values=B_component_array[dimension]
        )
        B_field[dimension] = interpolate_b(position_temp)

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


def get_grid_info(dimensions, D, flow_position, shock_index):
    if shock_index == 0:
        offset = 0

    # For plotting purposes arrays must have same dimensions.  Choose
    # plotting plane to have the same dimensions.
    if dimensions == 2:
        offset = 10
        if flow_position == 'upstream':
            if len(D.x2)/2 > shock_index:
                y_slice = len(D.x2)/2
            else:
                # don't take slice at shock position
                y_slice = shock_index + offset
        if flow_position == 'downstream':
            if len(D.x2)/2 < shock_index:
                y_slice = len(D.x2)/2
            else:
                # don't take slice at shock position
                y_slice = shock_index - offset
        grid = [D.x1, D.x3]
        B_component_array = [D.bx1[:, y_slice, :], D.bx3[:, y_slice, :]]
    if dimensions == 3:
        y_slice = -1  # not needed for 3D calculations
        # construct field lines a region that starts/ends a couple of grid
        # points beyond or before shock
        offset = 5
        if flow_position == 'upstream':
            grid = [D.x1, D.x2[shock_index + offset:], D.x3]
            B_component_array = [
                D.bx1[:, shock_index + offset:, :],
                D.bx2[:, shock_index + offset:, :],
                D.bx3[:, shock_index + offset:, :]
            ]
        if flow_position == 'downstream':
            grid = [D.x1, D.x2[:shock_index - offset + 1], D.x3]
            B_component_array = [
                D.bx1[:, :shock_index - offset + 1, :],
                D.bx2[:, :shock_index - offset + 1, :],
                D.bx3[:, :shock_index - offset + 1, :]
            ]

    nx = [0]*dimensions
    grid_min = [0.]*dimensions
    grid_max = [0.]*dimensions
    for dimension in range(0, dimensions):
        nx[dimension] = len(grid[dimension])
        grid_min[dimension] = min(grid[dimension])
        grid_max[dimension] = max(grid[dimension])

    return nx, grid, grid_min, grid_max, B_component_array, y_slice


def construct_field_lines(
    dimensions,
    D,
    separation,
    number_of_steps,
    step_length,
    flow_position,
    shock_index
):
    nx, grid, grid_min, grid_max, B_component_array, x_slice = \
        get_grid_info(dimensions, D, flow_position, shock_index)

    starting_positions = select_starting_positions(
        dimensions,
        nx,
        grid_min,
        grid_max,
        separation
    )
    # starting_positions = [[0.01, 0.124], [0.016, 0.124]]

    B_field = [[], []]
    for i in (0, 1):
        B_field[i] = interpolate_magnetic_field(
            dimensions,
            grid,
            B_component_array,
            starting_positions[i],
            shock_index
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
                grid_max,
                shock_index
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
                next_position,
                shock_index
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

    return [
        field_line_coordinates,
        field_line_components,
        step,
        x_slice,
    ]


def calculate_field_line_separation(
    dimensions,
    field_line_coordinates,
    field_line_components
):
    diffusion = [[], []]
    number_of_steps = len(field_line_coordinates[0][0])

    for step in range(0, number_of_steps):
        B_average = array([0.]*dimensions)
        separation_vector = array([0.]*dimensions)
        dot_product = 0.

        for dimension in range(0, dimensions):
            B_average[dimension] = (
                field_line_components[0][dimension][step] +
                field_line_components[1][dimension][step]
            )/2.

            separation_vector[dimension] = \
                field_line_coordinates[1][dimension][step] - \
                field_line_coordinates[0][dimension][step]

            dot_product = dot_product \
                + B_average[dimension]*separation_vector[dimension]

        B_average_magnitude = sqrt(sum(square(B_average)))   
        separation_magnitude = sqrt(sum(square(separation_vector)))
        theta = arccos(
            dot_product/(B_average_magnitude*separation_magnitude)
        )
        separation_perp_squared = square(separation_magnitude*sin(theta))

        diffusion[0].append(step)
        # diffusion[1].append(separation_magnitude_squared)
        diffusion[1].append(separation_perp_squared)

    return diffusion


def plot_field_lines(
    D,
    field_line_coordinates,
    field_line_components,
    x_slice,
    flow_position
):
    begin = 1
    end = 128

    field_line_coordinates_plot_y = field_line_coordinates[:][0]
    field_line_coordinates_plot_z = field_line_coordinates[:][1]
    for i in (0, 1):
        length = len(field_line_coordinates_plot_y[i])
        for j in range(0, length):
            y = field_line_coordinates_plot_y[i][j]
            z = field_line_coordinates_plot_z[i][j]
            if y < min(D.x2):
                y = max(D.x2) + (y - min(D.x2))
            if y > max(D.x2):
                y = min(D.x2) + (y - max(D.x2))
            if z < min(D.x3):
                z = max(D.x3) + (z - min(D.x3))
            if z > max(D.x3):
                z = min(D.x3) + (z - max(D.x3))
            field_line_coordinates_plot_y[i][j] = y
            field_line_coordinates_plot_z[i][j] = z

    axis('equal')
    quiver(
        D.x2[begin:end],
        D.x3[begin:end],
        D.bx2[x_slice, begin:end, begin:end].T,
        D.bx3[x_slice, begin:end, begin:end].T,
        units='xy'
    )

    colour = ['g', 'b']
    for i in (0, 0):
        axis('equal')
        quiver(
            field_line_coordinates_plot_y[0],
            field_line_coordinates_plot_z[1],
            field_line_components[i][0],
            field_line_components[i][1],
            units='xy',
            color=colour[i]
        )
        # axis([0, 0.025, 0, 0.025])
    show()


def add_result(number_of_steps, average_diffusion, diffusion):
    number_new_steps = len(diffusion[0])
    for step in range(0, number_new_steps):
        average_diffusion[1][step] = \
            average_diffusion[1][step] + diffusion[1][step]

    return average_diffusion


def write_diffusion_to_file(
    number_of_separations,
    average_diffusion_per_separation,
    flow_position,
    wdir
):
    f = open(wdir + flow_position + '_diffusion_field_lines.txt', 'w')
    # f.write('Distance along B,    RMS separation of lines\n')
    # f.write('\n')
    for separation in range(0, number_of_separations):
        number_of_steps = len(
            average_diffusion_per_separation[separation][0]
        )
        for step in range(0, number_of_steps):
            f.write(
                '%0.10f,        %1.10f\n' % (
                    average_diffusion_per_separation[separation][0][step],
                    average_diffusion_per_separation[separation][1][step]
                )
            )
        f.write('\n')  # NB to have another empty line


def calculate_B_field_line_diffusion(
    D,
    flow_position,
    shock_present,
    wdir
):
    dimensions = 3
    number_of_pairs = 1000
    number_of_steps = 1000  # along B
    initial_separations = [1e-3, 2e-3, 1e-2, 2e-2]
    number_of_separations = len(initial_separations)
    # assumes that the grid size is the same in all directions
    step_length = (max(D.x1) - min(D.x1))/(len(D.x1) - 1.)/1.
    max_plot_number = -1

    if shock_present:
        shock_index = locate_shock(D)
    else:
        shock_index = 0

    # set 'test = True' for testing
    test = False
    if test:
        dimensions = 2
        number_of_pairs = 1
        number_of_steps = 30  # along B
        step_length = (max(D.x1) - min(D.x1))/(len(D.x1) - 1.)/1.
        initial_separations = [1e-3]
        number_of_separations = len(initial_separations)
        max_plot_number = 4

    average_diffusion_per_separation = []
    for separation in initial_separations:
        number_usable_pairs = 0
        average_diffusion = [
            range(0, number_of_steps),
            [0]*number_of_steps
        ]
        for pair in range(0, number_of_pairs):
            print(separation, pair)
            field_line_coordinates, field_line_components, step, x_slice = \
                construct_field_lines(
                    dimensions,
                    D,
                    separation,
                    number_of_steps,
                    step_length,
                    flow_position,
                    shock_index
                )
            if step > -1:
                diffusion = calculate_field_line_separation(
                    dimensions,
                    field_line_coordinates,
                    field_line_components
                )
                add_result(number_of_steps, average_diffusion, diffusion)
                number_usable_pairs += 1

            # plot field lines - mainly for testing purposes
            plot_number = 1
            if dimensions == 2 and plot_number <= max_plot_number:
                plot_field_lines(
                    D,
                    field_line_coordinates,
                    field_line_components,
                    x_slice,
                    flow_position
                )
            plot_number += 1

        average_diffusion_per_separation.append([
            (array(average_diffusion[0]) + 1.)*step_length,
            sqrt(array(average_diffusion[1])/number_usable_pairs)
        ])

    write_diffusion_to_file(
        number_of_separations,
        average_diffusion_per_separation,
        flow_position,
        wdir
    )


if __name__ == '__main__':
    # flow_position: upstream or downstream of shock.  Options are:
    #                'upstream'
    #                'downstream'
    flow_position = 'upstream'
    file_number = 4
    wdir = '/home/mvorster/PLUTO/B_Shock_turbulence/bx0.8_cs1.35'
    shock_present = 1

    if not wdir[-1] == '/':
        wdir = wdir + '/'

    D = pload(file_number, w_dir=wdir+'output/')

    if not shock_present:
        flow_position = 'upstream'

    calculate_B_field_line_diffusion(
        D,
        flow_position,
        shock_present,
        wdir
    )
