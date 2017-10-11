# Turbulence analysis tools

A collection of tools to analyse the PLUTO simulation results.  The focus is specifically
on simulations where turbulence is present.


## Limitations
- Currently the code can only handle Cartesian coordinates    
- The code assumes that the simulations are done using the PLUTO fluid code
- Can only handle planar shock propagating along the x-axis


## Requirements
- Simulations run using PLUTO
- The `pyPLUTO.py` package that comes standard with PLUTO


## List of Tools

### `diffusion_field_lines.py`
- To investigate the diffusion of field lines in 3D:    
  * two random points in space are selected, separated by some distance 'l'    
  * the field lines are constructed in incremental steps starting from the two points     
  * the separation of the field lines, perpendicular to the local magnetic field, is calculated after every incremental step
  * the separation is calculated for a specified number of pairs of points     
- The results are written to a file    
- The diffusion can also be calculated in the presence of a planar shock    
  * assuming that the shock is propagating along the x-axis    
- For testing purposes field lines can also be plotted in 2D    

### `PLUTO_split_single_file.py`
- Splits a single PLUTO output file containing the data of all the fluid variables into separate files that 
only contain the data of a single variable

### `fitcurve.py`
- Calculates a fit to a curve using the Savitzky-Golan method

### `shock_tools.py`
- Contains a function to locate the position of the shock    
- Contains a function to calculate the turbulent variations (_delta_var_)    
  * the average fluid quantities are calculated along the propagation direction of the shock    
  * the average quantities are subtracted from the fluid fields to get the variations    

### `structure_function.py`
- Calculates the structure functions for the magnetic field and velocity
- Structure functions are calculated parallel and perpendicular to the local magnetic field
- Can handle the presence of a planar shock
  * can calculate the structure functions upstream and downstream of shock    
- Results are plotted    

### `units.py`
- Converts the program units used in PLUTO code to physical units

### `conserved_variables.py`
- compares the values of the flow variables (summed over the whole grid) at a specified simulation time with the initial values    
- for periodic boundary conditions the variables should be conserved

