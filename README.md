# Turbulence analysis tools

A collection of tools to analyse the PLUTO simulation results.  The focus is specifically
on simulations where turbulence is present, although the tools could be used for other
scnearios as well.


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
  * two random points in space are selected, separated by some distance _l_    
  * the field lines are constructed in incremental steps starting from the two points     
  * the separation of the field lines are calculated after every incremental step    
  * the separation is calculated for a specified number of pairs of points     
- The results are plotted and written to a file    
- The diffusion can alse be calculated in the presence of a planar shock    
  * assuming that the shock is propagating along the x-axis    
- For testing purposes field lines can also be plotted in 2D    
  * field lines are plotted    
  * it should in principle be easy to adapt code to construct field lines for other 2D simulations    

### `PLUTO_split_single_file.py`
- Splits a single PLUTO output file containing the data of all the fluid variables into separate files that 
only contain the data of a single variable

### `fitcurve.py`
- Calculates a fit to a curve using the Savitzky-Golan method

### `shock_tools.py`
- Contains a function to locate the position of the shock    
- Contains a function to calculate the turbulent variations (_delta_var_)    
  * the average fluid quantities along the x-direction are calculated    
  * (the shock is propagating along x-direction)    
  * the average quanties are subtracted from the fluid fields to get the variations    

### `structure_function.py`
- Calculates the structure functions for the density and velocity when a shock is present    
- Results are plotted    
- TODO: add for magnetic field    
- TODO remove unnecessary code   
- TODO calculate function upstream and downstream of shock

### `structure_function_initial.py`
- The same as previous module, but for the initial condition, i.e., where no shock is present

### `units.py`
- Converts the program units used in PLUTO code to physical units

