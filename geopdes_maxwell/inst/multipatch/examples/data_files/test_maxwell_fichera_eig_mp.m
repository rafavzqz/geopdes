% TEST_MAXWELL_FICHERA_EIG_MP: data file for Maxwell eigenproblem in Fichera's corner.

degree     = [2 2 2];     % Degree of the bsplines
regularity = [1 1 1];     % Regularity of the splines
n_sub      = [1 1 1];     % Number of subdivisions
nquad      = [3 3 3];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [];
drchlt_sides = [1];

% NURBS map from text file
geo_name = 'geo_fichera.txt';

% Physical parameters
c_elec_perm = @(x, y, z) ones(size(x));
c_magn_perm = @(x, y, z) ones(size(x));

