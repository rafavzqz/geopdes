% TEST_LINEAR_ELASTICITY_HORSESHOE: data file for linear elasticity problem on a horseshoe.

degree     = [3 3 3];     % Degree of the bsplines
regularity = [2 2 2];     % Regularity of the splines
n_sub      = [0 0 0];     % Number of subdivisions
nquad      = [4 4 4];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [1 2 3 4];
drchlt_sides = [5 6];
press_sides  = [];

% NURBS map from text file
geo_name = 'horseshoe.mat';

% Physical parameters
E  =  1; nu = 0; 
lam = @(x, y, z) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x))); 
mu  = @(x, y, z) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
fx = @(x, y, z) zeros (size (x));
fy = fx;
fz = @(x, y, z) ones (size (x));
f = @(x, y, z) cat(1, ...
                   reshape (fx (x,y,z), [1, size(x)]), ...
                   reshape (fy (x,y,z), [1, size(x)]), ...
                   reshape (fz (x,y,z), [1, size(x)]));
g       = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2));
h       = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2));

% Output file for Paraview
output_file = 'linear_elasticity_horseshoe';

% Points for post-processing
vtk_pts = {linspace(0, 1, 5)', linspace(0, 1, 5)', linspace(0, 1, 40)'};


% Exact solution
clear uex

