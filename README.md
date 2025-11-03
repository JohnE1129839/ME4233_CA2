# Overview
This repository is part of an assingment in Computational Methods in Fluid Dynamics. 
The code within this repository is meant to run a numerical analysis of the 2D flow within a container whose flow is subject to its lid velocity.
This simulation utilizes a numerical analysis on the stream-vorticity formulation of the Navier-Stokes equation to simulate the flow.
A Fourier analysis has also been conducted on a point of interest within the container to compare with the lid velocity frequency.

# Methods used
The following methods were used in the numerical analysis
- 5-point central difference scheme for spatial discretization
- Forward euler scheme for explicit time discretization
- Backward euler scheme for implicit time discretization

# How-to-use
This project does not require the use of outside add-ons and can run fully on MATLAB. To utilize this project follow the instructions below:
1. Set up the parameters in Parameters.txt
   - Re is the fluids Reynold's Nummber
   - Nx and Ny are the number of points in the x and y coordinates respectively
   - Lx and Ly are the 2-D dimensions of the container
   - dti is the timestep used for implicit numerical methods
   - dte is the timestep used for explicit numerical methods
   - tf is the simulation time duration
   - cx and cy are the position values for our point of interest
2. Set up the lid velocity profile in Unorth.m
3. Determine to use implicit or explicit methods in the lid_driven_cavity.m file
4. Run the lid_driven_cavity.m

# Results
The output of the code should include the following plots:
- A contour plot of the streamfunction at the final time
- A quiver plot of the flow velocity at the final time
- A plot about the horizontal flow velocity on the point of interest
- A Fourier analysis frequency domain plot on the point of interest
