% TEST_MAXWELL_CAMEMBERT: data file for Maxwell source problem in three quarters of the cylinder, where the exact solution is a singular function.

degree     = [2 2 2];     % Degree of the bsplines
regularity = [1 1 1];     % Regularity of the splines
n_sub      = [2 2 2];     % Number of subdivisions
nquad      = [3 3 3];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [3 4 5 6];
drchlt_sides = [1 2];

% NURBS map from text file
geo_name = 'geo_camembert.txt';

% Physical parameters
c_mass  = @(x, y, z) ones(size(x));
c_stiff = @(x, y, z) ones(size(x));

% Constant that characterizes the singularity
k = 1;

% Source and boundary terms
f       = @(x, y, z) cat (1, ...
                          singular_function (x, y, k), ...
                          zeros ([1, size(x)]));
g       = @(x, y, z, ind) zeros ([3, size(x)]);
h       = @(x, y, z, ind) zeros ([3, size(x)]);

% Exact solution
uex     = @(x, y, z) cat (1, ...
                          singular_function (x, y, k), ...
                          zeros ([1, size(x)]));
curluex = @(x, y, z) zeros (size(x));

% Output file for Paraview
output_file = 'maxwell_camembert';

% Points for post-processing
vtk_pts = {linspace(0, 1, 15)', linspace(0, 0.99, 15)', linspace(0, 1, 15)'};

% Remove warnings for the degenerated elements
warning ('off', 'geopdes:jacdet_zero_at_quad_node')
warning ('off', 'geopdes:zero_measure_element')

