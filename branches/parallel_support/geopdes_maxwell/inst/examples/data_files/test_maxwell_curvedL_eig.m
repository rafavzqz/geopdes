% TEST_MAXWELL_CURVEDL_EIG: data file for Maxwell eigenproblem in the curved L-shaped domain.

degree     = [3 3];     % Degree of the bsplines
regularity = [2 2];     % Regularity of the splines
n_sub      = [7 7];     % Number of subdivisions
nquad      = [4 4];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [];
drchlt_sides = [1 2 3 4];

% NURBS map from text file
geo_name = 'geo_curvedL.txt';

% Physical parameters
c_elec_perm = @(x, y) ones(size(x));
c_magn_perm = @(x, y) ones(size(x));

