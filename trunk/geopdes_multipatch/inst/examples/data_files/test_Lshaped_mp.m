% TEST_LSHAPED_MP: data file for Poisson problem in a multipatch L-shaped domain.

degree     = [3 3];  % Degree of the splines
regularity = [2 2];  % Regularity of the splines
n_sub      = [9 9];  % Number of subdivisions
nquad      = [5 5];  % Number of quadrature points

% NURBS map from a text file (in every case the result should be the same)
geo_name = 'geo_Lshaped_mp.txt';
%geo_name = 'geo_Lshaped_mp_b.txt';
%geo_name = 'geo_Lshaped_mp_c.txt';

% Type of boundary conditions.
nmnn_sides   = [3 4 5 6]; 
drchlt_sides = [1 2];

% Physical parameters
c_diff  = @(x, y) ones (size(x, 1), size(x, 2));

% Source and boundary terms
f = @(x, y) exp(x).*((x.^2 + y.^2 - 1).*sin(x.*y) - 2*y.*cos(x.*y));
g = @(x, y, ind) test_Lshaped_mp_g_nmnn (x, y, ind);
h = @(x, y, ind) exp(x) .* sin (x.*y);

% Exact solution
uex     = @(x, y) exp(x) .* sin (x.*y);
graduex = @(x, y) cat (1, ...
                       reshape ( exp(x).*(sin(x.*y) + y.*cos(x.*y)), [1, size(x)]), ...
                       reshape (exp(x).*x.*cos(x.*y), [1, size(x)]));

% Output file for Paraview
output_file = 'Lshaped_mp';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)'};
