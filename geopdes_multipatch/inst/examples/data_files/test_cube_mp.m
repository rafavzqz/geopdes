% TEST_CUBE_MP: data file for Poisson problem in a multipatch cube.

degree     = [2 2 2];  % Degree of the discrete functions
regularity = [1 1 1];  % Continuity of the discrete functions
n_sub      = [2 2 2];  % Number of knots for the mesh
nquad      = [3 3 3];  % Number of quadrature points

% Boundary conditions
drchlt_sides = [1];
nmnn_sides   = [];
if (~isempty (nmnn_sides))
  warning('test_cube_mp: we have not included the function test_g_nmnn for this case.')
end

% NURBS map from a text file
% You can see how the patches should match in the files
%  geo_2cubesa{b,c,d,e,f,g,h}.txt
% In all the eight cases the result must be the same
geo_name = 'geo_2cubesa.txt';

% Another geometry for testing
% geo_name = 'geo_4cubes.txt';

% Physical parameters
c_diff = @(x, y, z) ones (size(x));

% Source and boundary terms
f = @(x, y, z) -exp (x + z) .* sin (y);
h = @(x, y, z, ind) exp (x + z) .* sin (y); 
clear g

% Exact solution
uex     = @(x, y, z) exp (x + z) .* sin (y);
graduex = @(x, y, z) cat (1, ...
                          reshape (exp (x + z) .* sin (y), [1, size(x)]), ...
                          reshape (exp (x + z) .* cos (y), [1, size(x)]), ...
                          reshape (exp (x + z) .* sin (y), [1, size(x)]));

% Output file for Paraview
output_file = 'cube_mp';

% Points for post-processing
vtk_pts = {linspace(0, 1, 10)', linspace(0, 1, 15)', linspace(0, 1, 15)'};
