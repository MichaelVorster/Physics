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
# Last updated: 3 August 2018


import os
import pyPLUTO as pp
from fitcurve import savitzky_golay
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


def plot_shock_plane(D, shock_index_plane, shock_direction, wdir, file_time):
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    if shock_direction == 'x1':
        X = np.arange(D.x2[0], D.x2[-1]+D.dx2[0], D.dx2[0])
        Y = np.arange(D.x3[0], D.x3[-1]+D.dx3[0], D.dx3[0])
        Z = np.array(shock_index_plane)/float(D.n1)
        ax.set_xlabel('y')
        ax.set_ylabel('z')
        ax.set_zlabel('x')
        ax.set_zlim(0, D.x1r[-1])
    elif shock_direction == 'x2':
        X = np.arange(D.x1[0], D.x1[-1]+D.dx1[0], D.dx1[0])
        Y = np.arange(D.x3[0], D.x3[-1]+D.dx3[0], D.dx3[0])
        Z = np.array(shock_index_plane)/float(D.n2)
        ax.set_xlabel('x')
        ax.set_ylabel('z')
        ax.set_zlabel('y')
        ax.set_zlim(0, D.x2r[-1])
    else:
        X = np.arange(D.x1[0], D.x1[-1]+D.dx1[0], D.dx1[0])
        Y = np.arange(D.x2[0], D.x2[-1]+D.dx2[0], D.dx2[0])
        Z = np.array(shock_index_plane)/float(D.n3)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_zlim(0, D.x3r[-1])

    X, Y = np.meshgrid(X, Y)
    ax.plot_surface(X, Y, Z)
    plt.savefig(wdir + 'shock_plane_' + str(file_time) + '.png')
    # plt.show()
    plt.clf()


def locate_max_jump(pressure_slice, shock_pos_estimate):
    di = shock_pos_estimate/2
    i = shock_pos_estimate - di

    slice_avg = np.average(pressure_slice)

    while pressure_slice[i] > slice_avg:
        i += 1

    check_interval = range(i, i+di)
    for j in check_interval:
        if pressure_slice[j] > slice_avg:
            i = j + 1

    shock_index = i - 1

    return shock_index


def locate_shock_plane(
    D,
    shock_direction,
    shock_pos_estimate,
    wdir,
    file_time
):
    if shock_direction == 'x1':
        nx = [D.n2, D.n3]
    elif shock_direction == 'x2':
        nx = [D.n1, D.n3]
    else:
        nx = [D.n1, D.n2]

    shock_index_plane = np.zeros((nx[0], nx[1]), dtype=np.int_)

    for n0 in range(0, nx[0]):
        for n1 in range(0, nx[1]):
            if shock_direction == 'x1':
                pressure_slice = D.prs[:, n0, n1]
            elif shock_direction == 'x2':
                pressure_slice = D.prs[n0, :, n1]
            else:
                pressure_slice = D.prs[n0, n1, :]

            shock_index_plane[n0, n1] = locate_max_jump(
                pressure_slice,
                shock_pos_estimate
            )
            # shock_index_plane[n0, n1] = 246  # '+1' added few lines later

    plot_shock_plane(D, shock_index_plane, shock_direction, wdir, file_time)

    # Fortran is column-major order, and Python is row-major order.
    # np.array must therefore be transposed for Fortran purposes.
    # Should also add '1' to shock position since Fortran np.arrays start at 1
    shock_index_plane_transposed = np.transpose(shock_index_plane+1)
    output_file_name = wdir + 'shock_plane_' + str(file_time) + '.txt'
    np.savetxt(output_file_name, shock_index_plane_transposed)

    return np.average(shock_index_plane)


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

    plt.plot(x_grid, qx_slice, 'b', x_grid, avrg_qx, 'r', x_grid, avrg_fitted, 'k')
    plt.xlabel(shock_direction)
    plt.ylabel(r'Simulation [blue], np.average [red], Fitted [black]')
    plt.title('np.average ' + fluid_quantity)
    plt.axis([0, 1, y_min, y_max])
    plt.savefig(graph_dir + file_name + '_average')
    plt.show()


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
            avrg_qx[index] = np.average(qx[index, :, :])
        elif shock_direction == 'x2':
            avrg_qx[index] = np.average(qx[:, index, :])
        else:
            avrg_qx[index] = np.average(qx[:, :, index])

    avrg_fitted = list(avrg_qx)

    downstream_fit = savitzky_golay(
        avrg_fitted[0:shock_index], 41, 1
    )
    if shock_direction == 'x1':
        upstream_fit = np.average(qx[shock_index+1:, :, :])
    elif shock_direction == 'x2':
        upstream_fit = np.average(qx[:, shock_index+1:, :])
    else:
        upstream_fit = np.average(qx[:, :, shock_index+1:])

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
            del_qx.append(np.array(vx) - np.average(vx))
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
        del_qx.append(np.array(bx_para) - np.average(bx_para))
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

        del_component = np.array(list(qx[component]))
        for index in range(0, nx):
            if shock_direction == 'x1':
                del_component[index, :, :] = np.array(
                    qx[component][index, :, :] - avrg_qx[index]
                )
            elif shock_direction == 'x2':
                del_component[:, index, :] = np.array(
                    qx[component][:, index, :] - avrg_qx[index]
                )
            else:
                del_component[:, :, index] = np.array(
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
    wdir = '/home/mvorster/512_cube/bx1.6_cs0.8/1_perpendicular_shock/Run_1/PLUTO'
    file_time = 4
    shock_direction = 'x2'
    shock_pos_estimate = 256

    plot_averages = 0

    if not wdir[-1] == '/':
        wdir = wdir + '/'
    D = pp.pload(file_time, w_dir=wdir + 'output/')
    fluid_quantities = ['density', 'velocity', 'magnetic field']

    shock_index = locate_shock_plane(
      D,
      shock_direction,
      shock_pos_estimate,
      wdir,
      file_time
    )

    print(shock_index)

    # for fluid_quantity in fluid_quantities:
    #     remove_average_fluid_component(
    #         D,
    #         fluid_quantity,
    #         shock_index,
    #         shock_direction,
    #         wdir,
    #         plot_averages
    #     )
