# extrapolate_to_grid
Scripts for extrapolating discrete measurements such as temperature onto a hydrodynamic grid.

This code draws on `stompy` and provides a command line interface for extrapolating data
from a set of point measurements in time and space into a continuous in time and space
dataset. The spatial component of the extrapolation respects shorelines by performing the
extrapolation on the numerical grid.
