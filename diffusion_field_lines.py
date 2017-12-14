# Code to investigate diffusion of magnetic field lines.  Two field lines are
# constructed starting from two points separated by some initial distance.
# The separation between the field lines are calculated at set intervals along
# the lines.  This procedure is repeated for a number of field line pairs.
#
# AUTHOR: Michael Vorster
#
# LAST UPDATED: 14 December 2017


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
    shock_index,
    shock_direction
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

    # periodic boundaries, except in shock propagation direction
    
    if dimensions == 3:
        if shock_index > 0:
            if shock_direction == 'x1':
                i = 0
            if shock_direction == 'x2':
                i = 1
            if shock_direction == 'x3':
                i = 2
        if next_position[i] < grid_min[i] or next_position[i] > grid_max[i]:
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

    # periodic boundaries, except in shock propagation direction
    if shock_direction == 'x1':
        dimensions_range = [1, 2]
        position_temp[0] = position[0]
    if shock_direction == 'x2':
        dimensions_range = [0, 2]
        position_temp[1] = position[1]
    if shock_direction == 'x3':
        dimensions_range = [0, 1]
        position_temp[2] = position[2]

    if shock_index == 0 or dimensions == 2:
        dimensions_range = range(0, dimensions)

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


def get_grid_info(dimensions, D, flow_position, shock_index, shock_direction):
    if shock_index == 0:
        offset = 0

    # For plotting purposes arrays must have same dimensions.  Choose
    # plotting plane to have the same dimensions.
    if dimensions == 2:
        offset = 10
        len_coord = len(getattr(D, shock_direction))
        if flow_position == 'upstream':
            if len_coord/2 > shock_index + offset:
                grid_slice = len_coord/2
            else:
                # don't take slice at shock position
                grid_slice = shock_index + offset
        if flow_position == 'downstream':
            if len_coord/2 < shock_index - offset:
                grid_slice = len_coord/2
            else:
                # don't take slice at shock position
                grid_slice = shock_index - offset
        if shock_direction == 'x1':
            grid = [D.x2, D.x3]
            B_component_array = [
                D.bx2[grid_slice, :, :],
                D.bx3[grid_slice, :, :]
            ]
        elif shock_direction == 'x2':
            grid = [D.x1, D.x3]
            B_component_array = [
                D.bx1[:, grid_slice, :],
                D.bx3[:, grid_slice, :]
            ]
        else:
            grid = [D.x1, D.x2]
            B_component_array = [
                D.bx1[:, :, grid_slice],
                D.bx2[:, :, grid_slice]
            ]
    if dimensions == 3:
        grid_slice = -1  # not needed for 3D calculations
        # construct field lines a region that starts/ends a couple of grid
        # points beyond or before shock
        offset = 5
        if flow_position == 'upstream':
            cut_index = shock_index + offset
            xcut = cut_index if shock_direction == 'x1' else 0
            ycut = cut_index if shock_direction == 'x2' else 0
            zcut = cut_index if shock_direction == 'x3' else 0
            grid = [D.x1[xcut:], D.x2[ycut:], D.x3[zcut:]]
            B_component_array = [
                D.bx1[xcut:, ycut:, zcut:],
                D.bx2[xcut:, ycut:, zcut:],
                D.bx3[xcut:, ycut:, zcut:]
            ]
        if flow_position == 'downstream':
            cut_index = shock_index - offset + 1
            xcut = cut_index if shock_direction == 'x1' else D.n1
            ycut = cut_index if shock_direction == 'x2' else D.n2
            zcut = cut_index if shock_direction == 'x3' else D.n3
            grid = [D.x1[:xcut], D.x2[:ycut], D.x3[:zcut]]
            B_component_array = [
                D.bx1[:xcut, :ycut, :zcut],
                D.bx2[:xcut, :ycut, :zcut],
                D.bx3[:xcut, :ycut, :zcut]
            ]

    nx = [0]*dimensions
    grid_min = [0.]*dimensions
    grid_max = [0.]*dimensions
    for dimension in range(0, dimensions):
        nx[dimension] = len(grid[dimension])
        grid_min[dimension] = min(grid[dimension])
        grid_max[dimension] = max(grid[dimension])

    return nx, grid, grid_min, grid_max, B_component_array, grid_slice


def construct_field_lines(
    dimensions,
    D,
    separation,
    number_of_steps,
    step_length,
    flow_position,
    shock_index
):
    nx, grid, grid_min, grid_max, B_component_array, grid_slice = \
        get_grid_info(
            dimensions,
            D,
            flow_position,
            shock_index,
            shock_direction
        )

    starting_positions = select_starting_positions(
        dimensions,
        nx,
        grid_min,
        grid_max,
        separation
    )
    # starting_positions = [[0.1, 0.1], [0.12, 0.1]]

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
                shock_index,
                shock_direction
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
        grid_slice,
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
    grid_slice,
    flow_position,
    shock_direction
):
    begin = 45
    end = 70
    field_line_coordinates_plot_x = []
    field_line_coordinates_plot_y = []

    if shock_direction == 'x1':
        grid = ['x2', 'x3']
    if shock_direction == 'x2':
        grid = ['x1', 'x3']
    if shock_direction == 'x3':
        grid = ['x1', 'x2']

    for i in (0, 1):
        field_line_coordinates_plot_x.append(field_line_coordinates[i][0])
        field_line_coordinates_plot_y.append(field_line_coordinates[i][1])
        length = len(field_line_coordinates_plot_x[i])
        for j in range(0, length):
            x = field_line_coordinates_plot_x[i][j]
            y = field_line_coordinates_plot_y[i][j]
            xmin = min(getattr(D, grid[0]))
            xmax = max(getattr(D, grid[0]))
            ymin = min(getattr(D, grid[1]))
            ymax = max(getattr(D, grid[1]))

            if x < xmin:
                x = xmax + (x - xmin)
            if x > xmax:
                x = xmin + (x - xmax)
            if y < ymin:
                y = ymax + (y - ymin)
            if y > ymax:
                y = ymin + (y - ymax)
            field_line_coordinates_plot_x[i][j] = x
            field_line_coordinates_plot_y[i][j] = y

    axis('equal')
    if shock_direction == 'x1':
        quiver(
            getattr(D, grid[0])[begin:end],
            getattr(D, grid[1])[begin:end],
            D.bx2[grid_slice, begin:end, begin:end].T,
            D.bx3[grid_slice, begin:end, begin:end].T,
            units='xy'
        )
    elif shock_direction == 'x2':
        quiver(
            getattr(D, grid[0])[begin:end],
            getattr(D, grid[1])[begin:end],
            D.bx1[begin:end, grid_slice, begin:end].T,
            D.bx3[begin:end, grid_slice, begin:end].T,
            units='xy'
        )
    else:
        quiver(
            getattr(D, grid[0])[begin:end],
            getattr(D, grid[1])[begin:end],
            D.bx1[begin:end, begin:end, grid_slice].T,
            D.bx2[begin:end, begin:end, grid_slice].T,
            units='xy'
        )

    colour = ['g', 'b']
    for i in (0, 1):
        axis('equal')
        quiver(
            field_line_coordinates_plot_x[i],
            field_line_coordinates_plot_y[i],
            field_line_components[i][0],
            field_line_components[i][1],
            units='xy',
            color=colour[i]
        )
        #axis([0.1, 0.125, 0.115, 0.125])
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
    shock_direction,
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
        shock_index = locate_shock(D, shock_direction)
    else:
        shock_index = 0

    # set 'test = True' for testing
    test = False
    if test:
        dimensions = 2
        number_of_pairs = 1
        number_of_steps = 10  # along B
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
            field_line_coordinates, field_line_components, step, grid_slice = \
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
                    grid_slice,
                    flow_position,
                    shock_direction
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
    wdir = '/home/mvorster/512_cube/bx0.8_cs1.35/perpendicular_shock/Run_1/PLUTO'
    shock_present = 1
    shock_direction = 'x2'  # 'x1', 'x2', or 'x3'

    if not wdir[-1] == '/':
        wdir = wdir + '/'

    D = pload(file_number, w_dir=wdir+'output/')

    if not shock_present:
        flow_position = 'upstream'

    calculate_B_field_line_diffusion(
        D,
        flow_position,
        shock_present,
        shock_direction,
        wdir
    )
