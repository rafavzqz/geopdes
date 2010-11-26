% TEST_SQUARE_HEAT: data file for transient problem in the square.

degree     = [3 3];     % Degree of the bsplines
regularity = [2 2];     % Regularity of the splines
n_sub      = [9 9];     % Number of subdivisions
nquad      = [4 4];     % Points for the Gaussian quadrature rule

% Data for transient problem
deltat     = 1;
final_time = 20;
time_save  = linspace (0, final_time, final_time/deltat);

% Type of boundary conditions
nmnn_sides   = [];
drchlt_sides = [1 2 3 4];

% NURBS map from text file
geo_name = 'geo_square.txt';

% Physical parameters
c_mass  = @(x, y) ones(size(x));
c_diff  = @(x, y) ones(size(x));

% Source and boundary terms
f = @(x, y, t) 2*pi*pi*sin(pi*x).*sin(pi*y)*sin(t) + sin(pi*x).*sin(pi*y)*cos(t);
h = @(x, y, t, ind) sin(pi*x) .* sin(pi*y)*sin(t);

% Exact solution
uex     = @(x, y, t) sin(pi*x) .* sin(pi*y)*sin(t);
graduex = @(x, y, t) cat (1, ...
                     reshape (pi*cos(pi*x).*sin(pi*y)*sin(t), [1, size(x)]), ...
                     reshape (pi*sin(pi*x).*cos(pi*y)*sin(t), [1, size(x)]));

% Output file for Paraview
output_file = 'square_heat';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)'};

