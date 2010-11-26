% TEST_MAXWELL_CUBE: data file for Maxwell source problem in the cube.

degree     = [2 2 2];     % Degree of the bsplines
regularity = [1 1 1];     % Regularity of the splines
n_sub      = [2 2 2];     % Number of subdivisions
nquad      = [3 3 3];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [1 2];
drchlt_sides = [3 4 5 6];

% NURBS map from text file
geo_name = 'geo_cube.txt';

% Physical parameters
c_stiff = @(x, y, z) ones(size(x));
c_mass  = @(x, y, z) ones(size(x));

% Source and boundary terms
f       = @(x, y, z) cat(1, ...
                    reshape (sin(y) .* (2*z - exp(x)) - exp(z) .* sin(x), [1, size(x)]), ...
                    zeros ([1, size(x)]), ...
                    reshape (2 * exp(z) .* cos(x), [1, size(x)]));
g       = @(x, y, z, ind) test_maxwell_cube_g_nmnn (x, y, z, ind);
h       = @(x, y, z, ind) test_maxwell_cube_h_drchlt (x, y, z, ind);

% Exact solution
uex     = @(x, y, z) cat(1, ...
                    reshape (sin(y) .* z, [1, size(x)]), ...
                    reshape (exp(x) .* cos(y), [1, size(x)]), ...
                    reshape (exp(z) .* cos(x), [1, size(x)]));
curluex = @(x, y, z) cat(1, ...
                    zeros ([1, size(x)]), ...
                    reshape (sin(y) + exp(z) .* sin(x), [1, size(x)]), ...
                    reshape (exp(x).*cos(y) - z .* cos(y), [1, size(x)]));

% Output file for Paraview
output_file = 'maxwell_cube';

% Points for post-processing
vtk_pts = {linspace(0, 1, 15)', linspace(0, 1, 15)', linspace(0, 1, 15)'};
