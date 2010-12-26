% TEST_STOKES_3D_SYMDRIVCAV: data file for Stokes problem in the 3D driven cavity problem.

degree       = [2 2 2]; % Degree of the spline space for the pressure
regularity   = [1 1 1]; % Regularity of the spline space for the pressure
nquad        = [3 3 3]; % Number of quadrature points
nbreaks      = [3 3 3]; % Number of breaks (level of refinement)

% Description of the geometry
aspect_ratio = [1.5 .3];
geo_name     = vecscale ([1 aspect_ratio(1) aspect_ratio(2)]);

% Type of boundary conditions
drchlt_sides = 1:6;

% Discrete space (function to compute it)
der2         = false;

% Functions for the viscosity, the right-hand side, and boundary conditions
mu = @(x, y, z) ones (size (x));
f  = @(x, y, z) zeros ([3, size(x)]);
h  = @test_stokes_3d_symdrivcav_h_drchlt;

% Output file for Paraview
output_file  = 'test_stokes_3d_symdrivcav';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)', linspace(0, 1, 20)'};

clear uex graduex
