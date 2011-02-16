% TEST_PLANE_STRAIN_RING: data file for plane-strain problem on one quarter of a cylinder.

degree     = [3 3];     % Degree of the basis functions
regularity = [2 2];     % Regularity of the basis functions
n_sub      = [9 9];     % Number of subdivisions
nquad      = [4 4];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [2];
drchlt_sides = [];
press_sides  = [1];
symm_sides   = [3 4];

% NURBS map from text file
geo_name = 'geo_ring.txt';

% Physical parameters
E  =  1; nu = 0; 
lam = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x))); 
mu  = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
P = 1;
f = @(x, y) zeros (2, size (x, 1), size (x, 2));
g = @(x, y, ind) test_plane_strain_ring_g_nmnn (x, y, P, nu, ind);
h = @(x, y, ind) test_plane_strain_ring_uex (x, y, E, nu, P);
p = @(x, y, ind) P * ones (size (x));

% Exact solution
uex = @(x, y) test_plane_strain_ring_uex (x, y, E, nu, P);

% Output file for Paraview
output_file = 'plane_strain_ring';

% Points for post-processing
vtk_pts = {linspace(0, 1, 21)', linspace(0, 1, 21)'};

