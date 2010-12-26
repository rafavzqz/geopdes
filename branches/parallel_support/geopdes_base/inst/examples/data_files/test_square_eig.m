% TEST_SQUARE_EIG: data file for Laplace eigenproblem in the square.

degree     = [3 3];     % Degree of the bsplines
regularity = [2 2];     % Regularity of the splines
n_sub      = [7 7];     % Number of subdivisions
nquad      = [4 4];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [];
drchlt_sides = [1 2 3 4];

% NURBS map from text file
geo_name = 'geo_square.txt';

% Physical parameters
c_diff  = @(x, y) ones(size(x));
c_react = @(x, y) ones(size(x));
