% TEST_STOKES_SQUARE_BC_MP: data file for Stokes problem in the square defined with a multipatch geometry.

degree       = [3 3];   % Degree of the spline space for the pressure
regularity   = [2 2];   % Regularity of the spline space for the pressure
nquad        = [4 4];   % Number of quadrature points
nbreaks      = [10 10]; % Number of breaks (level of refinement)

% Type of boundary conditions
drchlt_sides = 1;

% Geometry description
geo_name = 'geo_square_multipatch3.txt';

% Function to compute the viscosity
mu = @(x, y) ones (size (x));

% Function to compute the right hand-side
fx = @(x, y) (6 * x + y .* cos(x .* y) + 2 * cos(y) .* sin(x));
fy = @(x, y) (x .* cos(x .* y) - 2 * cos(x) .* sin(y));

f = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));

% Exact solution, to compute the errors
uxex = @(x, y) (sin(x) .* cos(y));
uyex = @(x, y) (-sin(y) .* cos(x));

uex = @(x, y) cat(1, ...
                  reshape (uxex (x,y), [1, size(x)]), ...
                  reshape (uyex (x,y), [1, size(x)]));

pex = @(x, y) (3 * x.^2 + sin(x .* y) - 1.239811742000564725943866);

graduex = @test_stokes_square_bc_graduex;

% Functions to compute the Dirichlet boundary condition
h  = @(x, y, iside) uex (x, y);

% Output file for Paraview
output_file  = 'test_stokes_square_bc_mp';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)'};

