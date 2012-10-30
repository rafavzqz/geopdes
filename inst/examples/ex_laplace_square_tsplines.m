% EX_LAPLACE_SQUARE_TSPLINES: solve the Poisson problem in the unit square with a
%  T-spline discretization, coming from the T-splines plug-in for Rhino.

% PHYSICAL DATA OF THE PROBLEM
clear problem_data

% The physical domain is the unit square, defined either as a NURBS (square_structured) 
%  or as an unstructured T-spline.
% problem_data.geo_name = 'square_structured.iga'; % 2 x 2 structured mesh
problem_data.geo_name = 'square_unstructured.iga';

% Type of boundary conditions for each boundary set
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));

% Source and boundary terms
problem_data.f = @(x, y) 2 * pi^2 * sin(pi*x) .* sin (pi*y);
problem_data.h = @(x, y, ind) sin (pi*x) .* sin (pi*y);

% Exact solution (optional)
problem_data.uex     = @(x, y) sin (pi*x) .* sin (pi*y);
problem_data.graduex = @(x, y) cat (1, ...
                       reshape (pi*cos(pi*x).*sin(pi*y), [1, size(x)]), ...
                       reshape (pi*sin(pi*x).*cos(pi*y), [1, size(x)]));

% CALL TO THE SOLVER
[tspline, msh, space, u] = solve_laplace_tsplines_2d (problem_data);

% POST-PROCESSING
% Display errors of the computed solution in the L2 and H1 norm
[error_h1, error_l2] = ...
           tspline_h1_error (space, msh, u, problem_data.uex, problem_data.graduex)

% Export to Paraview
output_file = 'Square_tsplines';
vtk_pts_per_element = [10 10]; % Number of points on each Bezier element
tspline_to_vtk (u, tspline, vtk_pts_per_element, output_file, 'u')
