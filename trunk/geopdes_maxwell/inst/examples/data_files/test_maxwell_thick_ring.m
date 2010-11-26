% TEST_MAXWELL_THICK_RING: data file for Maxwell source problem in the thick ring.

degree     = [2 2 2];     % Degree of the bsplines
regularity = [1 1 1];     % Regularity of the splines
n_sub      = [2 2 2];     % Number of subdivisions
nquad      = [3 3 3];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [1 2];
drchlt_sides = [3 4 5 6];

% NURBS map from text file
geo_name = 'geo_thick_ring.txt';

% Physical parameters
c_stiff = @(x, y, z) ones(size(x));
c_mass  = @(x, y, z) ones(size(x));

% Source and boundary terms
f       = @(x, y, z) cat(1, ...
                    reshape (2*sin(y), [1, size(x)]), ...
                    reshape (2*sin(x), [1, size(x)]), ...
                    reshape (x .* y, [1, size(x)]));
g       = @(x, y, z, ind) test_maxwell_thick_ring_g_nmnn (x, y, z, ind);
h       = @(x, y, z, ind) test_maxwell_thick_ring_h_drchlt (x, y, z, ind);

% Exact solution
uex     = @(x, y, z) cat(1, ...
                    reshape (sin(y), [1, size(x)]), ...
                    reshape (sin(x), [1, size(x)]), ...
                    reshape (x .* y, [1, size(x)]));
curluex = @(x, y, z) cat(1, ...
                    reshape (x, [1, size(x)]), ...
                    reshape (-y, [1, size(x)]), ...
                    reshape (cos(x) - cos(y), [1, size(x)]));

% Output file for Paraview
output_file = 'maxwell_thick_ring';

% Points for post-processing
vtk_pts = {linspace(0, 1, 15)', linspace(0, 1, 15)', linspace(0, 1, 15)'};

