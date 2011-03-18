% TEST_MAXWELL_LSHAPED_EIG_MP: data file for Maxwell eigenproblem in the multipatch L-shaped domain.

degree     = [3 3];     % Degree of the bsplines
regularity = [2 2];     % Regularity of the splines
n_sub      = [5 5];     % Number of subdivisions
nquad      = [4 4];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [];
drchlt_sides = [1 2 3 4 5 6];

% NURBS map from text file (in the three cases the result should be the same)
geo_name = 'geo_Lshaped_mp.txt';
%geo_name = 'geo_Lshaped_mp_b.txt';
%geo_name = 'geo_Lshaped_mp_c.txt';

% Physical parameters
c_elec_perm = @(x, y) ones(size(x));
c_magn_perm = @(x, y) ones(size(x));

