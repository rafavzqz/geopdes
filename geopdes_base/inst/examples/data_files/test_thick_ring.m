% TEST_THICK_RING: data file for Poisson problem in a thick_ring.

degree     = [2 2 2];     % Degree of the bsplines
regularity = [1 1 1];     % Regularity of the splines
n_sub      = [3 3 3];     % Number of subdivisions
nquad      = [3 3 3];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [4 5 6];
drchlt_sides = [1 2 3];

% NURBS map from text file
geo_name = 'geo_thick_ring.txt';

% Physical parameters
c_diff  = @(x, y, z) ones(size(x));

% Source and boundary terms
f = @(x, y, z) exp(x).*cos(z).*(-2*cos(x.*y).*y + sin(x.*y).*(y.^2 + x.^2));
g = @test_thick_ring_g_nmnn;
h = @(x, y, z, ind) exp (x) .* sin (x.*y) .* cos (z);


% Exact solution
uex     = @(x, y, z) exp (x) .* sin (x.*y) .* cos (z);
graduex = @(x, y, z) cat (1, ...
                          reshape (exp (x) .* cos (z) .* (sin (x.*y) + y .* cos (x.*y)), [1, size(x)]), ...
                          reshape (exp (x) .* x .* cos (x.*y) .* cos(z), [1, size(x)]), ...
                          reshape (-exp (x) .* sin (x.*y) .* sin (z), [1, size(x)]));

% Output file for Paraview
output_file = 'thick_ring';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)', linspace(0, 1, 20)'};

