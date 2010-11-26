% TEST_STOKES_SQUARE: data file for Stokes problem in the square. The geometry is defined as an affine transformation.

degree       = [3 3];   % Degree of the spline space for the pressure
regularity   = [2 2];   % Regularity of the spline space for the pressure
nquad        = [4 4];   % Number of quadrature points
nbreaks      = [10 10]; % Number of breaks (level of refinement)

% Type of boundary conditions
drchlt_sides = 1:4;

% Geometry description
geo_name = eye(4);

% Discrete space (function to compute it)
sp_type = 'ndl';
switch sp_type
  case 'th'
    test_set_th_space;
  case 'ndl'
    test_set_ndl_space;
  case 'rt'
    test_set_rt_space;
end

% Functions to compute the viscosity and the Dirichlet boundary condition
mu = @(x, y) ones (size (x));
h  = @(x, y, iside) zeros ([2, size(x)]);

% Function to compute the right hand-side
fx = @(x, y) (exp(x).*(2.*(-1 + y).*y.*(2 + (-5 + y).*y) - ...
              8.*x.*(-1 + y).*y.*(2 + (-5 + y).*y) + ...
              6.*x.^3.*(-4 + (-5 + y).*(-2 + y).*y.*(1 + y)) + ...
              x.^2.*(12 + y.*(-38 + y.*(19 + (-6 + y).*y))) + ...
              x.^4.*(12 + y.*(-38 + y.*(19 + (-6 + y).*y)))));

fy = @(x, y) (456 - 912.*y + ...
              exp(x).*(-456 - 2.*x.*(-23 + 5.*x).*(10 + (-3 + x).*x) + ...
              2.*(456 + x.*(-466 + x.*(253 + x.*(-82 + 7.*x)))).* ...
              y  + (-6 + x.*(6 + x.*(-11 + x.*(22 + 7.*x)))).* ...
              y.^2 + 2.*(6 + x.*(10 + x.*(-29 + (-6 + x).*x))).* ...
              y.^3 + (-6 + x.*(3 + x).*(-2 + x.*(7 + x))).*y.^4));

f = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));

% Exact solution, to compute the errors
uxex = @(x, y) (2.*exp(x).*(-1 + x).^2.*x.^2.*(-1 + y).*y.*(-1 + 2.*y));
uyex = @(x, y) (-exp(x).*(-1 + x).*x.*(-2 + x.*(3 + x)).*(-1 + y).^2.*y.^2);

uex = @(x, y) cat(1, ...
                  reshape (uxex (x,y), [1, size(x)]), ...
                  reshape (uyex (x,y), [1, size(x)]));

pex = @(x, y) (-424 + 156.*exp(1) + (-1 + y).*y.* ...
              (-456 + exp(x).*(456 + x.^2.*(228 - 5.*(-1 + y).*y) ...
               + 2.*x.*(-228 + (-1 + y).*y) ...
               + 2.*x.^3.*(-36 + (-1 + y).*y) ...
               + x.^4.*(12 + (-1 + y).*y))));

graduex = @test_stokes_square_graduex;

% Output file for Paraview
output_file  = 'test_stokes_square';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)'};

