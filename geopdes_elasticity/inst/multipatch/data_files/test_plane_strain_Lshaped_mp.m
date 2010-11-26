% TEST_PLANE_STRAIN_LSHAPED_MP: data file for plane strain problem on a multipatch L-shaped domain.

degree     = [3 3];     % Degree of the bsplines
regularity = [2 2];     % Regularity of the splines
n_sub      = [5 5];     % Number of subdivisions
nquad      = [4 4];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [];
drchlt_sides = [1 2 3 4 5 6];

% NURBS map from text file
% In the three cases the result must be the same
geo_name = 'geo_Lshaped_mp.txt';
% geo_name = 'geo_Lshaped_mp_b.txt';
% geo_name = 'geo_Lshaped_mp_c.txt';

% Physical parameters
E  =  1; nu = 0.3;
lam = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x))); 
mu  = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
fx = @(x, y) lam (x, y) .* sin(y) + mu(x, y) .* (1 + x) .* sin(y);
fy = @(x, y) mu (x, y) .* cos(y) .* (2*x - 1) - lam (x, y) .* cos(y) .* (1 - x);
f  = @(x, y) cat(1, ...
                   reshape (fx (x,y), [1, size(x)]), ...
                   reshape (fy (x,y), [1, size(x)]));
h = @(x, y, ind) cat(1, ...
                        reshape (x.*sin(y), [1, size(x)]), ...
			reshape (x.*cos(y), [1, size(x)]));
clear g p

% Exact solution
uxex = @(x, y) x.*sin(y);
uyex = @(x, y) x.*cos(y);
uex  = @(x, y) cat(1, ...
		      reshape (uxex (x,y), [1, size(x)]), ...
		      reshape (uyex (x,y), [1, size(x)]));

% Output file for Paraview
output_file = 'plane_strain_Lshaped_mp';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)'};

