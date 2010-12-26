% TEST_MAXWELL_THICKL_EIG: data file for Maxwell eigenproblem in the thick L-shaped domain.

degree     = [2 2 2];     % Degree of the bsplines
regularity = [1 1 1];     % Regularity of the splines
n_sub      = [2 2 2];     % Number of subdivisions
nquad      = [3 3 3];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [];
drchlt_sides = [1 2 3 4 5 6];

% NURBS map from text file
geo_name = 'geo_thickL_C0.txt';
%geo_name = 'geo_thickL_C1.txt';

% Physical parameters
c_elec_perm = @(x, y, z) ones(size(x));
c_magn_perm = @(x, y, z) ones(size(x));

