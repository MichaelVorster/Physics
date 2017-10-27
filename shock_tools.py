# Calculate the average profiles of fluid quantities along the propagation axis
# of the planar shock.  The fluid quantities are averaged in the plane perpen-
# dicular to the propagation direction.  In general the averaged profile will
# not be smooth.  Therefore, it is necessary to fit a smooth function to the
# profile.  However, the presence of a shock complicates the problem.
#
# The approach is to fit the shocked region with a smooth curve.  Upstream of
# the shock the fluid quantity is again averaged, but this time over all axes.
# The second averaging should ensure a smooth profile upstream of the shock.

# Author: Michael Vorster
# Last updated: 10 October 2017


import os
import pyPLUTO as pp
from fitcurve import savitzky_golay
from numpy import (
    argmax,
    array,
    average,
    int_,
    square,
    sqrt,
    zeros
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


def create_scalar_fields(D):
    density = D.rho
    velocity = sqrt(
        square(D.vx1) +
        square(D.vx2) +
        square(D.vx3)
    )
    magnetic_field = sqrt(
        square(D.bx1) +
        square(D.bx2) +
        square(D.bx3)
    )

    pressure = D.prs
    # scalar_fields = [
    #     density,
    #     velocity,
    #     magnetic_field
    # ]

    # using just pressure for bx0.8_cs1.35 seems sufficient
    scalar_fields = [
        pressure
    ]

    return scalar_fields


def locate_max_jump(variables):
    avrg_index_variables = []
    for variable in variables:
        n = 10
        dn = 0
        dn_max = 100
        initial_estimate = argmax(variable)
        ratio = []
        while dn < dn_max:
            set_1 = variable[
                initial_estimate - n + dn: initial_estimate + dn + 1
            ]
            set_2 = variable[
                initial_estimate + dn + 1: initial_estimate + n + dn
            ]
            ratio.append(average(set_1)/average(set_2))
            dn += 1
        avrg_index_variables.append(
            initial_estimate + argmax(ratio)
        )
    avrg_index = int(average(avrg_index_variables))

    return avrg_index


def locate_shock(D, shock_direction):
    variables = create_scalar_fields(D)
    nx = len(getattr(D, shock_direction))

    # average over plane perpendicular to shock direction
    avrg_variables = []
    for variable in variables:
        avrg_qx = [0]*nx
        for index in range(0, nx):
            if shock_direction == 'x1':
                avrg_qx[index] = average(variable[index, :, :])
            elif shock_direction == 'x2':
                avrg_qx[index] = average(variable[:, index, :])
            else:
                avrg_qx[index] = average(variable[:, :, index])
        avrg_variables.append(avrg_qx)

    shock_index = locate_max_jump(avrg_variables)

    return shock_index


def locate_shock_plane(D, shock_direction):
    variables = create_scalar_fields(D)

    if shock_direction == 'x1':
        nx = [D.n2, D.n3]
    elif shock_direction == 'x2':
        nx = [D.n1, D.n3]
    else:
        nx = [D.n1, D.n2]

    shock_index_plane = zeros((nx[0], nx[1]), dtype=int_)

    for n0 in range(0, nx[0]):
        for n1 in range(0, nx[1]):
            variables_slices = []
            for variable in variables:
                if shock_direction == 'x1':
                    single_variable_slice = variable[:, n0, n1]
                elif shock_direction == 'x2':
                    single_variable_slice = variable[n0, :, n1]
                else:
                    single_variable_slice = variable[n0, n1, :]
                variables_slices.append(single_variable_slice)

            shock_index_plane[n0, n1] = locate_max_jump(variables_slices)

    return shock_index_plane


def plot_average_fluid_field(
    D,
    x_grid,
    qx,
    avrg_qx,
    avrg_fitted,
    fluid_quantity,
    shock_direction,
    wdir
):
    graph_dir = wdir + 'fit_to_components_averaged/'
    if not os.path.exists(graph_dir):
        os.makedirs(graph_dir)

    if shock_direction == 'x1':
        qx_slice = qx[:, len(D.x2)/2, len(D.x3)/2]
    elif shock_direction == 'x2':
        qx_slice = qx[len(D.x1)/2, :, len(D.x3)/2]
    else:
        qx_slice = qx[len(D.x1)/2, len(D.x2)/2, :]

    y_min = 1.2*min(qx_slice)
    y_max = 1.2*max(qx_slice)
    file_name = fluid_quantity[:]
    for char in ['$', '{', '}']:
        file_name = file_name.replace(char, '')

    plot(x_grid, qx_slice, 'b', x_grid, avrg_qx, 'r', x_grid, avrg_fitted, 'k')
    xlabel(shock_direction)
    ylabel(r'Simulation [blue], Average [red], Fitted [black]')
    title('Average ' + fluid_quantity)
    axis([0, 1, y_min, y_max])
    savefig(graph_dir + file_name + '_average')
    show()


def average_fluid_field(
    D,
    x_grid,
    qx,
    nx,
    fluid_quantity,
    shock_index,
    shock_direction,
    wdir,
    plot_averages
):
    avrg_qx = [0]*nx
    for index in range(0, nx):
        if shock_direction == 'x1':
            avrg_qx[index] = average(qx[index, :, :])
        elif shock_direction == 'x2':
            avrg_qx[index] = average(qx[:, index, :])
        else:
            avrg_qx[index] = average(qx[:, :, index])

    avrg_fitted = list(avrg_qx)

    downstream_fit = savitzky_golay(
        avrg_fitted[0:shock_index], 41, 1
    )
    if shock_direction == 'x1':
        upstream_fit = average(qx[shock_index+1:, :, :])
    elif shock_direction == 'x2':
        upstream_fit = average(qx[:, shock_index+1:, :])
    else:
        upstream_fit = average(qx[:, :, shock_index+1:])

    for index in range(0, shock_index):
        avrg_fitted[index] = downstream_fit[index]
    for index in range(shock_index, nx):
        avrg_fitted[index] = upstream_fit

    if plot_averages:
        plot_average_fluid_field(
            D,
            x_grid,
            qx,
            avrg_qx,
            avrg_fitted,
            fluid_quantity,
            shock_direction,
            wdir
        )

    return avrg_fitted


def remove_average_fluid_component(
    D,
    fluid_quantity,
    shock_index,
    shock_direction,
    wdir,
    plot_averages
):
    nx = len(getattr(D, shock_direction))
    del_qx = []

    if shock_direction == 'x1':
        perp_components = ['x2', 'x3']
    elif shock_direction == 'x2':
        perp_components = ['x1', 'x3']
    else:
        perp_components = ['x1', 'x2']

    if fluid_quantity == 'density':
        qx = [D.rho]
        qx_title = ['density']
    elif fluid_quantity == 'velocity':
        qx = [getattr(D, 'v' + shock_direction)]
        qx_title = ['$V_{' + shock_direction + '}$']
        for perp_component in perp_components:
            vx = getattr(D, 'v' + perp_component)
            del_qx.append(array(vx) - average(vx))
    elif fluid_quantity == 'magnetic field':
        qx = [
            getattr(D, 'b' + perp_components[0]),
            getattr(D, 'b' + perp_components[1])
        ]
        qx_title = [
            '$B_{' + perp_components[0] + '}$',
            '$B_{' + perp_components[1] + '}$'
        ]
        bx_para = getattr(D, 'b' + shock_direction)
        del_qx.append(array(bx_para) - average(bx_para))
    else:
        print('Fluid quantity not recognised')
        exit(0)

    number_of_components = len(qx)
    for component in range(0, number_of_components):
        avrg_qx = average_fluid_field(
            D,
            getattr(D, shock_direction),
            qx[component],
            nx,
            qx_title[component],
            shock_index,
            shock_direction,
            wdir,
            plot_averages
        )

        del_component = array(list(qx[component]))
        for index in range(0, nx):
            if shock_direction == 'x1':
                del_component[index, :, :] = array(
                    qx[component][index, :, :] - avrg_qx[index]
                )
            elif shock_direction == 'x2':
                del_component[:, index, :] = array(
                    qx[component][:, index, :] - avrg_qx[index]
                )
            else:
                del_component[:, :, index] = array(
                    qx[component][:, :, index] - avrg_qx[index]
                )
        del_qx.append(del_component)

    if fluid_quantity == 'velocity':
        if shock_direction == 'x1':
            del_qx[0], del_qx[1], del_qx[2] = del_qx[2], del_qx[0], del_qx[1]
        if shock_direction == 'x2':
            del_qx[1], del_qx[2] = del_qx[2], del_qx[1]

    return del_qx


if __name__ == '__main__':
    # wdir = '/home/mvorster/PLUTO/Shock_turbulence/Results/Run_15_b/'
    # file_time = 7
    # shock_direction = 'x1'
    # shock_index = 535
    wdir = '/home/mvorster/512_cube/bx0.8_cs1.35/Run_1/PLUTO/'
    file_time = 4
    shock_direction = 'x2'
    shock_index = 260

    plot_averages = 1

    if not wdir[-1] == '/':
        wdir = wdir + '/'
    D = pp.pload(file_time, w_dir=wdir + 'output/')
    fluid_quantities = ['density', 'velocity', 'magnetic field']

    shock_index = locate_shock_plane(D, shock_direction)
    # print(shock_index)
    # for fluid_quantity in fluid_quantities:
    #     remove_average_fluid_component(
    #         D,
    #         fluid_quantity,
    #         shock_index,
    #         shock_direction,
    #         wdir,
    #         plot_averages
    #     )
