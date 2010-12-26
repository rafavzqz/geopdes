% TEST_RING: data file for Poisson problem in one quarter of a ring, with homogeneous Dirichlet boundary conditions.

degree     = [3 3];     % Degree of the bsplines
regularity = [2 2];     % Regularity of the splines
n_sub      = [9 9];     % Number of new knots on each interval
nquad      = [4 4];     % Points for the Gaussian quadrature rule

% Type of boundary conditions
nmnn_sides   = [];
drchlt_sides = [1 2 3 4];

% NURBS map from text file (.txt) or a Matlab binary file (.mat). Both should give the same result
geo_name = 'geo_ring.txt';
% geo_name = 'ring.mat'

% Physical parameters
c_diff  = @(x, y) ones(size(x));

% Source and boundary terms
f       = @(x, y) 2*x.*(22.*x.^2.*y.^2+21.*y.^4-45.*y.^2+x.^4-5.*x.^2+4);
h       = @(x, y, ind) zeros (size (x));
clear g

% Exact solution
uex     = @(x, y) -(x.^2+y.^2-1).*(x.^2+y.^2.-4).*x.*y.^2;
graduex = @(x, y) cat (1, ...
                       reshape (-2*(x.*y).^2.*((x.^2+y.^2-1)+(x.^2+y.^2-4)) - ...
                                (x.^2+y.^2-1).*(x.^2+y.^2-4).*y.^2, [1, size(x)]), ...
                       reshape ( -2*x.*y.^3.*((x.^2+y.^2-1)+(x.^2+y.^2-4)) - ...
                                2*x.*y.*(x.^2+y.^2-1).*(x.^2+y.^2-4), [1, size(x)]));

% Output file for Paraview
output_file = 'ring';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)'};

