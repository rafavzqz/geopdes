% TEST_MAXWELL_PACMAN: data file for Maxwell source problem in three quarters of a circle, where the exact solution is a singular function.

degree     = [3 3];     % Degree of the bsplines
regularity = [2 2];     % Regularity of the splines
n_sub      = [5 5];     % Number of subdivisions
nquad      = [4 4];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [3 4];
drchlt_sides = [1 2];

% NURBS map from text file
geo_name = 'geo_pacman.txt';


% Physical parameters
c_mass  = @(x, y) ones(size(x));
c_stiff = @(x, y) ones(size(x));

% Constant that characterizes the singularity
k = 1;

% Source and boundary terms
f       = @(x, y) singular_function (x, y, k);
g       = @(x, y, ind) zeros([2, size(x)]);
h       = @(x, y, ind) zeros (size(x));

% Exact solution
uex     = @(x, y) singular_function (x, y, k);
curluex = @(x, y) zeros (size(x));

% Output file for Paraview
output_file = 'maxwell_pacman';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, .99, 20)'};

% Remove warnings for the degenerated elements
warning ('off', 'geopdes:jacdet_zero_at_quad_node')
warning ('off', 'geopdes:zero_measure_element')
