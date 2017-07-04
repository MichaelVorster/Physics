# Simulations consist of a rectangular grid, with x-axis along the long-axis
# of rectangular region.  A planar shock propagates along x-axis
#
# Calculate the average profile of a fluid quantity along the x-axis, i.e.,
# the fluid quantity is averaged along the y- and z-axes.  In general the
# averaged profile will not be smooth.  Therefore, it is necessary to fit a
# smooth function to the profile.  However, the presence of a shock complicates
# the problem.
#
# The approach is to fit the shocked region with a smooth curve.  Upstream of
# the shock the fluid quantity is again averaged, but this time over all axes.
# The second averaging should ensure a smooth profile upstream of the shock.


import pyPLUTO as pp
from fitcurve import savitzky_golay
from numpy import (
    argmax,
    array,
    average,
    square,
    sqrt
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


def locate_shock(D):
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
    variables = [
        density,
        velocity,
        magnetic_field
    ]
    nx = len(D.x1)

    # average over y and z directions
    avrg_variables = []
    for variable in variables:
        avrg_qx = [0]*nx
        for i in range(0, nx):
            avrg_qx[i] = average(variable[i, :, :])
        avrg_variables.append(avrg_qx)

    avrg_shock_index_variables = []
    for avrg_variable in avrg_variables:
        shock_index = []
        n = 3
        dn = 0
        dn_max = 30
        initial_estimate = argmax(avrg_variable)
        ratio = []
        while dn < dn_max:
            set_1 = avrg_variable[
                initial_estimate - n + dn : initial_estimate + dn + 1  # noqa
            ]
            set_2 = avrg_variable[
                initial_estimate + dn + 1 : initial_estimate + n + dn  # noqa
            ]
            ratio.append(average(set_1)/average(set_2))
            dn += 1
        shock_index.append(initial_estimate + argmax(ratio))
        avrg_shock_index_variables.append(int(average(shock_index)) + 1)
    avrg_shock_index = int(average(avrg_shock_index_variables)) + 1

    return avrg_shock_index


def grid_variables(D):
    nx = (len(D.x1), len(D.x2), len(D.x3))
    dx = [
        (max(D.x1) - min(D.x1))/(nx[0] - 1.),
        (max(D.x2) - min(D.x2))/(nx[1] - 1.),
        (max(D.x3) - min(D.x3))/(nx[2] - 1.)
    ]

    return nx, dx


def plot_average_fluid_field(x1, qx, avrg_qx, avrg_fitted, fluid_quantity):
    y_max = 1.2*max(avrg_qx)
    plot(x1, qx[:, 80, 80], 'b', x1, avrg_qx, 'r', x1, avrg_fitted, 'k')
    xlabel(r'x')
    ylabel(r'Simulation [blue], Average [red], Fitted [black]')
    title('Average ' + fluid_quantity)
    axis([0, 1, -0.5, y_max])
    savefig(fluid_quantity+'_x_average')
    show()


def average_fluid_field(x1, qx, nx, fluid_quantity, shock_index):
    avrg_qx = [0]*nx[0]
    for i in range(0, nx[0]):
        avrg_qx[i] = average(qx[i, :, :])

    avrg_fitted = list(avrg_qx)

    downstream_fit = savitzky_golay(
        avrg_fitted[0:shock_index], 41, 1
    )
    upstream_fit = average(qx[shock_index+1:, :, :])

    for i in range(0, shock_index):
        avrg_fitted[i] = downstream_fit[i]
    for i in range(shock_index, nx[0]):
        avrg_fitted[i] = upstream_fit

    plot_average_fluid_field(x1, qx, avrg_qx, avrg_fitted, fluid_quantity)

    return avrg_fitted


def remove_average_fluid_component(D, fluid_quantity, shock_index):
    nx, dx = grid_variables(D)
    del_qx = []

    if flow_quantity == 'density':
        qx = [D.rho]
        qx_title = ['density']
    elif flow_quantity == 'velocity':
        qx = [D.vx1]
        qx_title = ['$V_x$']
        del_qx.append(array(D.vx2) - average(D.vx2))
        del_qx.append(array(D.vx3) - average(D.vx3))
    elif flow_quantity == 'magnetic field':
        qx = [D.bx2, D.bx3]
        qx_title = ['$B_y$', '$B_z$']
        del_qx.append(array(D.bx1) - average(D.bx1))
    else:
        print('Fluid quantity not recognised')
        exit(0)

    number_of_components = len(qx)
    for component in range(0, number_of_components):
        avrg_qx = average_fluid_field(
            D.x1,
            qx[component],
            nx,
            qx_title[component],
            shock_index
        )

        del_component = array(list(qx[component]))
        for i in range(0, nx[0]):
            del_component[i, :, :] = array(qx[component][i, :, :] - avrg_qx[i])
        del_qx.append(del_component)

    if flow_quantity == 'velocity':
        del_qx[0], del_qx[1], del_qx[2] = del_qx[2], del_qx[0], del_qx[1]

    return del_qx


if __name__ == '__main__':
    wdir = '/home/mvorster/PLUTO/Shock_turbulence/Results/Run_15/output/'
    file_time = 9
    D = pp.pload(file_time, w_dir=wdir)

    # flow_quantity: options are
    #                'density'
    #                'velocity'
    #                'magnetic field'
    flow_quantity = 'magnetic field'

    # shock_index = locate_shock(D)
    shock_index = 502
    remove_average_fluid_component(D, flow_quantity, shock_index)
