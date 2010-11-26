% TEST_THICK_LSHAPED_MP: data file for Poisson problem in a multipatch thick L-shaped domain.

degree     = [2 2 2];  % Degree of the discrete functions
regularity = [1 1 1];  % Continuity of the discrete functions
n_sub      = [2 2 2];  % Number of knots for the mesh
nquad      = [3 3 3];  % Number of quadrature points

% Nurbs map from a text file (in both cases the result should be the same)
geo_name = 'geo_thickL_mp.txt';
%geo_name = 'geo_thickL_mp_b.txt';

% Type of boundary conditions
nmnn_sides   = [1 2]; 
drchlt_sides = [3 4 5 6 7 8];

% Physical parameters
c_diff = @(x, y, z) ones (size(x));

% Source and boundary terms
f = @(x, y, z) exp(x).*cos(z).*(-2*cos(x.*y).*y + sin(x.*y).*(y.^2+x.^2));
g = @(x, y, z, ind) test_thick_Lshaped_mp_g_nmnn (x, y, z, ind);
h = @(x, y, z, ind) exp (x) .* sin (x.*y) .* cos (z);

% Exact solution
uex     = @(x, y, z) exp (x) .* sin (x.*y) .* cos (z);
graduex = @(x, y, z) cat (1, ...
                          reshape (exp (x) .* cos (z) .* (sin (x.*y) + y .* cos (x.*y)), [1, size(x)]), ...
                          reshape (exp (x) .* x .* cos (x.*y) .* cos(z), [1, size(x)]), ...
                          reshape (-exp (x) .* sin (x.*y) .* sin (z), [1, size(x)]));

% Output file for Paraview
output_file = 'test_thick_Lshaped_mp';

% Points for post-processing
vtk_pts = {linspace(0, 1, 10)', linspace(0, 1, 10)', linspace(0, 1, 20)'};
