% TEST_SQUARE_HEAT_EXP: data file for transient problem in the square.

degree     = [3 3];     % Degree of the bsplines
regularity = [2 2];     % Regularity of the splines
n_sub      = [9 9];     % Number of subdivisions
nquad      = [4 4];     % Points for the Gaussian quadrature rule

% Data for transient problem
deltat     = 1;
final_time = 5;
time_save  = linspace (0, final_time, final_time/deltat);

% Type of boundary conditions
nmnn_sides   = [1 2 3];
drchlt_sides = [4];

% NURBS map from text file
geo_name = 'geo_square.txt';

% Physical parameters
c_mass  = @(x, y) ones(size(x));
c_diff  = @(x, y) ones(size(x));

% Source and boundary terms
f = @(x, y, t) -2 * exp(y) * exp(t);
g = @(x, y, t, ind) test_square_heat_exp_g_nmnn (x, y, t, ind);
h = @(x, y, t, ind) x .* x .* exp(y) * exp(t);

% Exact solution
uex     = @(x, y, t) x .* x .* exp(y) * exp(t);
graduex = @(x, y, t) cat (1, ...
                     reshape (2*x.*exp(y)*exp(t), [1, size(x)]), ...
                     reshape (x.*x.*exp(y)*exp(t), [1, size(x)]));

% Output file for Paraview
output_file = 'square_heat';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)'};

