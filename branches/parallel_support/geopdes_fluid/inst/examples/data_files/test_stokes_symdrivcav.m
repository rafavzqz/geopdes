% TEST_STOKES_SYMDRIVCAV: data file for Stokes problem in the driven cavity problem.

degree       = [2 2];   % Degree of the spline space for the pressure
regularity   = [1 1];   % Regularity of the spline space for the pressure
nquad        = [4 4];   % Number of quadrature points
nbreaks      = [10 10]; % Number of breaks (level of refinement)

% Description of the geometry
aspect_ratio = 0.9;
geo_name     = vecscale ([1 aspect_ratio 1]);

% Type of boundary conditions
drchlt_sides = 1:4;

% Discrete space (function to compute it)
sp_type = 'th';
switch sp_type
  case 'th'
    test_set_th_space; 
  case 'ndl'
    test_set_ndl_space;
  case 'rt'
    test_set_rt_space;
end

% Functions for the viscosity, the right-hand side, and boundary conditions
mu = @(x, y) ones (size (x));
f  = @(x, y) zeros ([2, size(x)]);
h  = @test_stokes_symdrivcav_h_drchlt;

% Output file for Paraview
output_file  = 'test_stokes_symdrivcav';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)'};


clear uex graduex
