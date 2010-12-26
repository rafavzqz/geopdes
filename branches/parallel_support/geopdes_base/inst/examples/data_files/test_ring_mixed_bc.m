% TEST_RING_MIXED_BC: data file for Poisson problem in one quarter of a ring.

degree     = [3 3];     % Degree of the bsplines
regularity = [2 2];     % Regularity of the splines
n_sub      = [8 8];     % Number of subdivisions
nquad      = [4 4];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [1 3 4];
drchlt_sides = [2];

% NURBS map from text file
geo_name = 'geo_ring.txt';

% Physical parameters
c_diff  = @(x, y) ones(size(x));

% Source and boundary terms
f = @(x, y) exp(x).*((x.^2 + y.^2 - 1).*sin(x.*y) - 2*y.*cos(x.*y));
g = @test_ring_mixed_bc_g_nmnn;
h = @(x, y, ind) exp(x).*sin(x.*y);


% Exact solution
uex     = @(x, y) exp(x) .* sin (x.*y);
graduex = @(x, y) cat (1, ...
                       reshape ( exp(x).*(sin(x.*y) + y.*cos(x.*y)), [1, size(x)]), ...
                       reshape (exp(x).*x.*cos(x.*y), [1, size(x)]));

% Output file for Paraview
output_file = 'ring_mixed_bc';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)'};

