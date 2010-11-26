% TEST_CUBE: data file for Poisson problem in the unit cube.

degree     = [3 3 3];     % Degree of the bsplines
regularity = [2 2 2];     % Regularity of the splines
n_sub      = [2 2 2];     % Number of new knots on each interval
nquad      = [4 4 4];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [4 5 6];
drchlt_sides = [1 2 3];

% NURBS map from text file
geo_name = 'geo_cube.txt';

% Physical parameters
c_diff  = @(x, y, z) ones(size(x));

% Source and boundary terms
f = @(x, y, z) -exp(x+z).*sin(y);
g = @test_cube_g_nmnn;
h = @(x, y, z, ind) exp (x + z) .* sin (y);

% Exact solution
uex     = @(x, y, z) exp (x + z) .* sin (y);
graduex = @(x, y, z) cat (1, ...
                          reshape (exp (x + z) .* sin (y), [1, size(x)]), ...
                          reshape (exp (x + z) .* cos (y), [1, size(x)]), ...
                          reshape (exp (x + z) .* sin (y), [1, size(x)]));

% Output file for Paraview
output_file = 'cube';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)', linspace(0, 1, 20)'};

