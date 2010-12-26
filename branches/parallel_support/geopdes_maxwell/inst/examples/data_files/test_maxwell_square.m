% TEST_MAXWELL_SQUARE: data file for Maxwell source problem in the square.

degree     = [3 3];     % Degree of the bsplines
regularity = [2 2];     % Regularity of the splines
n_sub      = [7 7];     % Number of subdivisions
nquad      = [4 4];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [2 3];
drchlt_sides = [1 4];

% NURBS map from text file
geo_name = 'geo_square.txt';

% Physical parameters
c_stiff = @(x, y) ones(size(x));
c_mass  = @(x, y) ones(size(x));

% Source and boundary terms
f       = @(x, y) cat(1, ...
                    reshape (-exp(x) .* sin(y) + 2*sin(y), [1, size(x)]), ...
                    zeros ([1, size(x)]));
g       = @(x, y, ind) test_maxwell_square_g_nmnn (x, y, ind);
h       = @(x, y, ind) test_maxwell_square_h_drchlt (x, y, ind);

% Exact solution
uex     = @(x, y) cat(1, ...
                    reshape (sin(y), [1, size(x)]), ...
                    reshape (exp(x) .* cos(y), [1, size(x)]));
curluex = @(x, y) exp(x).*cos(y) - cos(y);

% Output file for Paraview
output_file = 'maxwell_square';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)'};

