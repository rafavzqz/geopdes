% TEST_PLANE_STRAIN_PLATE: data file for plane-strain problem on a square plate with a circular hole.

degree     = [3 3];     % Degree of the basis functions
regularity = [2 2];     % Regularity of the basis functions
n_sub      = [7 7];     % Number of subdivisions
nquad      = [4 4];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [];
drchlt_sides = [3];
press_sides  = [1 2 4];
symm_sides   = [];

% NURBS map from text file
geo_name = 'geo_plate_with_hole.txt';

% Physical parameters
E  =  1; nu = .3; 
lam = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x))); 
mu  = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
f = @(x, y) zeros (2, size (x, 1), size (x, 2));
h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
p = @(x, y, ind) ones (size (x));

% Output file for Paraview
output_file = 'plane_strain_plate';

% Points for post-processing
vtk_pts = {linspace(0, 1, 21)', linspace(0, 1, 21)'};


clear uex g
