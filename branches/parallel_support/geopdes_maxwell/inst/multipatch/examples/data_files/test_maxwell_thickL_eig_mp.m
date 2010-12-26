% TEST_MAXWELL_THICKL_EIG_MP: data file for Maxwell eigenproblem in the multipatch thick L-shaped domain.

degree     = [2 2 2];     % Degree of the bsplines
regularity = [1 1 1];     % Regularity of the splines
n_sub      = [2 2 2];     % Number of subdivisions
nquad      = [3 3 3];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [];
drchlt_sides = [1 2 3 4 5 6 7 8];

% NURBS map from text file (in both cases the result should be the same)
geo_name = 'geo_thickL_mp.txt';
%geo_name = 'geo_thickL_mp_b.txt';

% Physical parameters
c_elec_perm = @(x, y, z) ones(size(x));
c_magn_perm = @(x, y, z) ones(size(x));

