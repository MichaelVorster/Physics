# Code to investigate diffusion of magnetic field lines.  Two field lines are
# constructed starting from two points separated by some initial distance.
# The separation between the field lines are calculated at set intervals along
# the lines.  This procedure is repeated for a number of field line pairs.
#
# AUTHOR: Michael Vorster
#
# LAST UPDATED: 07 September 2017


from math import (
    sin,
    cos
)
from matplotlib.pyplot import (
    axis,
    loglog,
    quiver,
    savefig,
    show,
    xlabel,
    ylabel
)
from numpy import (
    array,
    mod,
    power,
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
    file_number
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

    # periodic boundaries in y and z, and in x when t = 0
    #   (just be careful about the pressure wave at inner
    #   boundary for PLUTO simulations)
    if file_number > 0:
        if next_position[0] < grid_min[0] or next_position[0] > grid_max[0]:
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
    file_number
):
    B_field = [0.]*dimensions
    position_temp = [0.]*dimensions

    # periodic boundaries along y and z (and x at t = 0) are taken into account
    dimensions_range = range(1, dimensions)
    if file_number == 0 or dimensions == 2:
        dimensions_range = range(0, dimensions)

    # periodic boundaries along y and z (and x at t = 0) are taken into account
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


def get_grid_info(file_number, dimensions, D, flow_position, shock_index):
    if file_number == 0:
        offset = 0
        if flow_position == 'upstream':
            shock_index = 0
        if flow_position == 'downstream':
            shock_index = max(D.x1)

    # for plotting purposes arrays must have same dimensions.
    # Only x2 and x3 have the same dimensions. Slice is thus taken in
    # x-y plane.
    if dimensions == 2:
        offset = 10
        if flow_position == 'upstream':
            if len(D.x1)/2 > shock_index:
                x_slice = len(D.x1)/2
            else:
                # don't take slice at shock position
                x_slice = shock_index + offset
        if flow_position == 'downstream':
            if len(D.x1)/2 < shock_index:
                x_slice = len(D.x1)/2
            else:
                # don't take slice at shock position
                x_slice = shock_index - offset
        grid = [D.x2, D.x3]
        B_component_array = [D.bx2[x_slice, :, :], D.bx3[x_slice, :, :]]
    if dimensions == 3:
        x_slice = -1  # not needed for 3D calculations
        # construct field lines a region that starts/ends a couple of grid
        # points beyond or before shock
        offset = 5
        if flow_position == 'upstream':
            grid = [D.x1[shock_index + offset:], D.x2, D.x3]
            B_component_array = [
                D.bx1[shock_index + offset:, :, :],
                D.bx2[shock_index + offset:, :, :],
                D.bx3[shock_index + offset:, :, :]
            ]
        if flow_position == 'downstream':
            grid = [D.x1[:shock_index - offset + 1], D.x2, D.x3]
            B_component_array = [
                D.bx1[:shock_index - offset + 1, :, :],
                D.bx2[:shock_index - offset + 1, :, :],
                D.bx3[:shock_index - offset + 1, :, :]
            ]

    nx = [0]*dimensions
    grid_min = [0.]*dimensions
    grid_max = [0.]*dimensions
    for dimension in range(0, dimensions):
        nx[dimension] = len(grid[dimension])
        grid_min[dimension] = min(grid[dimension])
        grid_max[dimension] = max(grid[dimension])

    return nx, grid, grid_min, grid_max, B_component_array, x_slice


def construct_field_lines(
    file_number,
    dimensions,
    D,
    separation,
    number_of_steps,
    step_length,
    flow_position,
    shock_index
):
    nx, grid, grid_min, grid_max, B_component_array, x_slice = \
        get_grid_info(file_number, dimensions, D, flow_position, shock_index)

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
            file_number
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
                file_number
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
                file_number
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


def calculate_field_line_separation(dimensions, field_lines):
    diffusion = [[], []]
    number_of_steps = len(field_lines[0][0])
    for step in range(0, number_of_steps):
        separation_vector = array([0.]*dimensions)
        for dimension in range(0, dimensions):
            separation_vector[dimension] = \
                field_lines[1][dimension][step] - \
                field_lines[0][dimension][step]
        separation_magnitude_squared = sum(square(separation_vector))
        '''
        if dimensions == 2:
            separation_perp_squared = square(separation_vector[1])
        if dimensions == 3:
            separation_perp_squared = (
                square(separation_vector[0]) +
                square(separation_vector[2])
            )
        '''
        diffusion[0].append(step)
        diffusion[1].append(separation_magnitude_squared)
        # diffusion[1].append(separation_perp_squared)

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


def plot_diffusion(
    number_of_separations,
    average_diffusion_per_separation,
    flow_position
):
    Richardson_x = array(range(10, 100))*1e-3
    Richardson_y = power(Richardson_x, 1.5)*0.5
    Kolmogorov_x = array(range(10, 100))*1e-3
    Kolmogorov_y = power(Kolmogorov_x, 0.5)*0.2
    for separation in range(0, number_of_separations):
        loglog(
            average_diffusion_per_separation[separation][0],
            average_diffusion_per_separation[separation][1]
        )
    loglog(Richardson_x, Richardson_y)
    loglog(Kolmogorov_x, Kolmogorov_y)
    xlabel(r'Distance along B')
    ylabel(r'RMS separation of lines')
    savefig(flow_position + '_B_field_line_diffusion.png')
    # show()


def write_diffusion_to_file(
    number_of_separations,
    average_diffusion_per_separation,
    flow_position
):
    f = open(flow_position + '_diffusion_field_lines.txt', 'w')
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


def calculate_B_field_line_diffusion(file_number, D, flow_position):
    dimensions = 3
    number_of_pairs = 1000
    number_of_steps = 1000  # along B
    initial_separations = [1e-4, 5e-4, 1e-3, 5e-3]
    number_of_separations = len(initial_separations)
    # assumes that the grid size is the same in all directions
    step_length = (max(D.x1) - min(D.x1))/(len(D.x1) - 1.)/2.
    max_plot_number = -1

    if file_number > 0:
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
                    file_number,
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
                    field_line_coordinates
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
        flow_position
    )
    plot_diffusion(
        number_of_separations,
        average_diffusion_per_separation,
        flow_position
    )


if __name__ == '__main__':
    # flow_position: upstream or downstream of shock.  Options are:
    #                'upstream'
    #                'downstream'
    flow_position = 'downstream'
    file_number = 0
    # wdir = '/home/mvorster/PLUTO/Shock_turbulence/Results/Run_15_b/output/'
    wdir = '/home/cronus/vorster/PLUTO/Shock_turbulence/turbulence_only_output/'

    D = pload(file_number, w_dir=wdir)

    if file_number == 0:
        flow_position = 'upstream'  # no shock at time = 0

    calculate_B_field_line_diffusion(file_number, D, flow_position)
