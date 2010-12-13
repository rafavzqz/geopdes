% TEST_WAVEGUIDE_3D_MP: data file for linear elasticity problem on a multipatch waveguide.

degree     = [2 2 2];     % Degree of the bsplines
regularity = [1 1 1];     % Regularity of the splines
n_sub      = [1 1 1];     % Number of subdivisions
nquad      = [3 3 3];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [];
drchlt_sides = [3 4];

% NURBS map from text file
% You can see how the patches should match in the files
%  geo_2cubesa{b,c,d,e,f,g,h}.txt
% In all the eight cases the result must be the same
geo_name = 'geo_square_section_waveguide_mp.txt';

% Physical parameters
E  =  1; nu = 0.3;
lam = @(x, y, z) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x))); 
mu  = @(x, y, z) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
fx = @(x, y, z) 0*x;
fy = @(x, y, z) 0*x;
fz = @(x, y, z) 0*x;
f = @(x, y, z) cat(1, ...
                   reshape (fx (x,y,z), [1, size(x)]), ...
                   reshape (fy (x,y,z), [1, size(x)]), ...
                   reshape (fz (x,y,z), [1, size(x)]));
h = @test_waveguide_3d_mp_h_drchlt;
clear g p uex

% Output file for Paraview
output_file = 'linear_elasticity_waveguide_mp';

% Points for post-processing
vtk_pts = {linspace(0, 1, 10)', linspace(0, 1, 10)', linspace(0, 1, 10)'};

