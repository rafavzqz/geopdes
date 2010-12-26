% TEST_PLATE_MIXED_BC: data file for Poisson problem in one quarter of the plate with hole (but without symmetry).

degree     = [3 3];     % Degree of the bsplines
regularity = [2 2];     % Regularity of the splines
n_sub      = [7 7];   % Number of new knots on each interval
nquad      = [4 4];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [3 4];
drchlt_sides = [1 2];

% NURBS map from text file
geo_name = 'geo_plate_with_hole.txt';

% Physical parameters
c_diff  = @(x, y) ones(size(x));

% Source and boundary terms
f = @(x, y) zeros(size(x));
g = @test_plate_mixed_bc_g_nmnn;
h = @(x, y, ind)  exp(x).*sin(y);


% Exact solution
uex     = @(x, y) exp(x).*sin (y);
graduex = @(x, y) cat (1, ...
                       reshape (exp(x).*sin(y), [1, size(x)]), ...
                       reshape (exp(x).*cos(y), [1, size(x)]));

% Output file for Paraview
output_file = 'plate_mixed_bc';

% Points for post-processing
vtk_pts = {linspace(0, 1, 21)', linspace(0, 1, 21)'};

