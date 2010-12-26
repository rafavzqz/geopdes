% TEST_MAXWELL_CUBE_EIG: data file for Maxwell eigenproblem in the cube.

degree     = [2 2 2];     % Degree of the bsplines
regularity = [1 1 1];     % Regularity of the splines
n_sub      = [3 3 3];     % Number of subdivisions
nquad      = [3 3 3];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [];
drchlt_sides = [1 2 3 4 5 6];

% NURBS map from text file
geo_name = 'geo_cube.txt';

% Physical parameters
c_elec_perm = @(x, y, z) ones(size(x));
c_magn_perm = @(x, y, z) ones(size(x));

