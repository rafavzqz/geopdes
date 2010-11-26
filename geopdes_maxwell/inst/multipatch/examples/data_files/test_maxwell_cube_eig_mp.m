% TEST_MAXWELL_CUBE_EIG_MP: data file for Maxwell eigenproblem in a multipatch cube.

degree     = [2 2 2];     % Degree of the bsplines
regularity = [1 1 1];     % Regularity of the splines
n_sub      = [1 1 1];     % Number of subdivisions
nquad      = [3 3 3];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [];
drchlt_sides = [1];

% NURBS map from text file
% You can see how the patches should match in the files 
%  geo_2cubesa{b,c,d,e,f,g,h}.txt
% In all the eight cases the result must be the same
geo_name = 'geo_2cubesa.txt';

% Another geometry for testing
% geo_name = 'geo_4cubes.txt';

% Physical parameters
c_elec_perm = @(x, y, z) ones(size(x));
c_magn_perm = @(x, y, z) ones(size(x));

