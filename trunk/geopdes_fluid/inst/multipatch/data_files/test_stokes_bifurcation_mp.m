% TEST_STOKES_BIFURCATION_MP: data file for Stokes problem in the square defined with a multipatch geometry.

degree       = [3 3];   % Degree of the spline space for the pressure
regularity   = [2 2];   % Regularity of the spline space for the pressure
nquad        = [4 4];   % Number of quadrature points
nbreaks      = [10 10]; % Number of breaks (level of refinement)

% Type of boundary conditions
drchlt_sides = 1:3;

% Geometry description
geo_name = 'geo_bifurcation_mp.txt';

% Function to compute the viscosity
mu = @(x, y) ones (size (x));

% Function to compute the right hand-side
fx = @(x, y) zeros (size (x));
fy = @(x, y) zeros (size (x));

f = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));

% Functions to compute the Dirichlet boundary condition
h  = @(x, y, iside) uex (x, y);

% Output file for Paraview
output_file  = 'test_stokes_bifurcation_mp';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)'};

