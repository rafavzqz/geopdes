% TEST_PLANE_STRAIN_SQUARE: data file for plane-strain problem on a square.

degree     = [3 3];     % Degree of the bsplines
regularity = [2 2];     % Regularity of the splines
n_sub      = [9 9];     % Number of subdivisions
nquad      = [4 4];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [];
press_sides  = [];
drchlt_sides = [1 2 3 4];
symm_sides   = [];

% NURBS map constructed using the NURBS toolbox
geo_name = nrb4surf([0 0], [1 0], [0 1], [1 1]);

% Physical parameters
E  =  1; nu = .3; 
lam = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x))); 
mu  = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
fx = @(x, y) -(-(mu(x, y)*3 + lam(x, y)).*sin(2*pi*x).*sin(2*pi*y) + (mu(x, y) + lam(x, y)).*cos(2*pi*x).*cos(2*pi*y))*(2*pi)^2;
fy = fx;
f = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));
h       = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
clear g p

% Exact solution
uxex = @(x,y) sin(2*pi*x).*(sin(2*pi*y));
uyex = @(x,y) sin(2*pi*x).*(sin(2*pi*y));
uex = @(x, y) cat(1, ...
                reshape (uxex (x,y), [1, size(x)]), ...
                reshape (uyex (x,y), [1, size(x)]));

% Output filename for Paraview
output_file = 'plane_strain_square';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)'};

