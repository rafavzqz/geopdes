% TEST_SQUARE: data file for Poisson problem in the square.

degree     = [3 3];     % Degree of the bsplines
regularity = [2 2];     % Regularity of the splines
n_sub      = [9 9];     % Number of subdivisions
nquad      = [4 4];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [];
drchlt_sides = [1 2 3 4];

% NURBS map from text file
geo_name = 'geo_square.txt';

% Physical parameters
c_diff  = @(x, y) ones(size(x));

% Source and boundary terms
f = @(x, y) zeros (size (x));
g = @test_square_g_nmnn;
h = @(x, y, ind) exp (x) .* sin(y);

% Exact solution
uex     = @(x, y) exp (x) .* sin (y);
graduex = @(x, y) cat (1, ...
                       reshape (exp(x).*sin(y), [1, size(x)]), ...
                       reshape (exp(x).*cos(y), [1, size(x)]));

% Output file for Paraview
output_file = 'square';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)'};

