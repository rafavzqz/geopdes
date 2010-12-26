% TEST_MAXWELL_RING: data file for Maxwell source problem in one quarter of a ring.

degree     = [3 3];     % Degree of the bsplines
regularity = [2 2];     % Regularity of the splines
n_sub      = [9 9];     % Number of subdivisions
nquad      = [4 4];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [];
drchlt_sides = [1 2 3 4];

% NURBS map from text file
geo_name = 'geo_ring.txt';

% Physical parameters
c_stiff = @(x, y) ones(size(x));
c_mass  = @(x, y) ones(size(x));

% Source and boundary terms
f       = @(x, y) cat(1, ...
                    reshape (2*sin(y), [1, size(x)]), ...
                    reshape (2*sin(x), [1, size(x)]));
g       = @(x, y, ind) test_maxwell_ring_g_nmnn (x, y, ind);
h       = @(x, y, ind) test_maxwell_ring_h_drchlt (x, y, ind);

% Exact solution
uex     = @(x, y) cat(1, ...
                    reshape (sin(y), [1, size(x)]), ...
                    reshape (sin(x), [1, size(x)]));
curluex = @(x, y) cos(x) - cos(y);

% Output file for Paraview
output_file = 'maxwell_ring';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)'};

