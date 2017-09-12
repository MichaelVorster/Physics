# This code is used to split a single PLUTO file containing the data for all
# variables into a number of single files that contain the data for only one
# variable


import array
import numpy as np


data_dir = '/home/mvorster/PLUTO/Shock_turbulence/Results/Run_22'
file_time = '0000'

input_file_name = 'data.' + file_time + '.dbl'
input_file_path = data_dir + '/output/' + input_file_name

endian = "<"
dtype = 'd'
n1 = 1024
n2 = 128
n3 = 128

input_file = open(input_file_path, 'rb')
myvars = [
    'rho',
    'vx1',
    'vx2',
    'vx3',
    'bx1',
    'bx2',
    'bx3',
    'prs',
    'psi_glm'
]

fmt = endian+str(n1*n2*n3)+dtype
nb = np.dtype(fmt).itemsize

for myvar in myvars:
    A = array.array(dtype)
    A.fromstring(input_file.read(nb))

    output_file_name = myvar + '.' + file_time + '.dbl'
    output_file_path = data_dir + '/particle_code_input/' + output_file_name
    output_file = open(output_file_path, 'wb')
    A.tofile(output_file)
