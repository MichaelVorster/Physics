# The primary purpose of the code is to determine how the conserved fluid
# variables (density, momentum, magnetic field, and energy) vary with time.
# The results at some simulation time is compared to the initial values.
# For a fully periodic system the variables should be conserved.


import collections
from numpy import square, sqrt, sum
from Tools import pyPLUTO as pp


wdir = '/home/mvorster/PLUTO/Shock_turbulence/Results/Run_21/output/'
file_times = [0, 1]

# equation of state: 'adiabatic' or 'isothermal'
EOS = 'adiabatic'

# B_field_present: 1: = MHD or 0 = HD
B_field_present = 1


if EOS == 'adiabatic':
    # adiabatic index
    gamma = 5./3.
if EOS == 'isothermal':
    sound_speed = 5./3.

rho = []
prs = []
V_x = []
V_y = []
V_z = []
B_x = []
B_y = []
B_z = []

V = []
B = []
p_x = []
p_y = []
p_z = []
kinetic_energy = []
internal_energy = []
magnetic_energy = []
energy = []


i = 0
for file_time in file_times:
    D = pp.pload(file_time, w_dir=wdir)

    rho.append(D.rho)
    if EOS == 'adiabatic':
        prs.append(D.prs)
    if EOS == 'isothermal':
        prs.append(square(sound_speed)*D.rho)
    V_x.append(D.vx1)
    V_y.append(D.vx2)
    V_z.append(D.vx3)
    if B_field_present:
        B_x.append(D.bx1)
        B_y.append(D.bx2)
        B_z.append(D.bx3)
    else:
        B_x.append(1e-12)
        B_y.append(1e-12)
        B_z.append(1e-12)

    V.append(
        sqrt(
            square(V_x[i]) +
            square(V_y[i]) +
            square(V_z[i])
        )
    )
    B.append(
        sqrt(
            square(B_x[i]) +
            square(B_y[i]) +
            square(B_z[i])
        )
    )
    p_x.append(D.rho*D.vx1)
    p_y.append(D.rho*D.vx2)
    p_z.append(D.rho*D.vx3)
    kinetic_energy.append(0.5*D.rho*V[i]*V[i])
    if EOS == 'adiabatic':
        internal_energy.append(prs[i]/(gamma - 1.))
    if EOS == 'isothermal':
        internal_energy.append(prs[i]/prs[i]*1e-12)
    magnetic_energy.append(0.5*B[i]*B[i])
    energy.append(kinetic_energy[i] + internal_energy[i] + magnetic_energy[i])

    i += 1

variables = collections.OrderedDict([
    ('Density:', rho),
    ('p_x:', p_x),
    ('p_y:', p_y),
    ('p_z:', p_z),
    ('B_x:', B_x),
    ('B_y:', B_y),
    ('B_z:', B_z),
    ('Energy:', energy),

    ('V:', V),
    ('V_x:', V_x),
    ('V_y:', V_y),
    ('V_z:', V_z),
    ('B:', B),
    ('Pressure:', prs),
    ('Kinetic energy:', kinetic_energy),
    ('Magnetic energy:', magnetic_energy),
    ('Internal energy:', internal_energy)
])

# write to file
f = open(wdir + 'conserved_variables.txt', 'w')
format_string = '{0:17} {1:1.12f}\n'

for file_time in file_times:
    f.write('t = ' + str(file_time) + '\n')
    f.write('=====\n')

    for variable in variables:
        f.write(
            format_string.format(
                variable,
                sum(variables[variable][file_time])
            )
        )
        if variable == 'Energy:':
            f.write('\n')

    f.write('\n')
    f.write('\n')

f.write('Relative difference\n')
f.write('===================\n')

for variable in variables:
    init_val = sum(variables[variable][0])
    final_val = sum(variables[variable][1])
    rel_diff = (final_val - init_val)/init_val
    f.write(
        format_string.format(
            variable,
            rel_diff
        )
    )
    if variable == 'Energy:':
        f.write('\n')

f.close()
