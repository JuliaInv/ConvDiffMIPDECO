# Finite Difference Code ConvDiffMIPDECO

This folder contains AMPL models of MIPDECO instances for the two-dimensional source inversion problem. For simplicity we use a finite difference discretization of the PDE. 

## Overview

1. `SrcInvertFDM.mod` is the main ampl model
1. `solver.ampl` selects the solver
1. `runFile.sh` is a script that runs the experiment
1. `SrcInvertFDMN10.dat` contains the measurements
1. The files `N08.ampl`, `N16.ampl`, and `N32.ampl` are the models for grid sizes of 8x8, 16x16, and 32x32, respectively
1. `plotUW.ampl` creates a MATLAB script to plot the solution

