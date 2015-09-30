% EX_LAPLACE_SQUARE: solve the Poisson problem in the unit square with a B-spline discretization.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = nrbline ([0 0], [1 0]);

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2];

% Physical parameters
problem_data.c_diff  = @(x) ones(size(x));

% Source and boundary terms
n = 10;
problem_data.f = @(x) n^2 * pi^2 * sin (n*pi*x);
problem_data.g = @(x, ind) n * pi * cos (n*pi*x) * (2*ind-3); % The last part fixes the sign
problem_data.h = @(x, ind) sin (n*pi*x);

% Exact solution (optional)
problem_data.uex     = @(x) sin (n*pi*x);
problem_data.graduex = @(x) n * pi * cos (n*pi*x);

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = 3;     % Degree of the splines
method_data.regularity = 2;     % Regularity of the splines
method_data.nsub       = 20;    % Number of subdivisions
method_data.nquad      = 4;     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_laplace (problem_data, method_data);

% 4) POST-PROCESSING
% PLOT IN MATLAB. COMPARISON WITH THE EXACT SOLUTION

vtk_pts = {linspace(0, 1, 150)};
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
plot (F, eu, F, problem_data.uex(F))
legend ('Numerical solution', 'Exact solution'); axis tight

% Display errors of the computed solution in the L2 and H1 norm
[error_h1, error_l2] = ...
           sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex)

%!demo
%! ex_laplace_1d

%!test
%! problem_data.geo_name = nrbline ([0 0], [1 0]);
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2];
%! problem_data.c_diff  = @(x) ones(size(x));
%! n = 10;
%! problem_data.f = @(x) n^2 * pi^2 * sin (n*pi*x);
%! problem_data.g = @(x, ind) n * pi * cos (n*pi*x) * (2*ind-3); % The last part fixes the sign
%! problem_data.h = @(x, ind) sin (n*pi*x);
%! problem_data.uex     = @(x) sin (n*pi*x);
%! problem_data.graduex = @(x) n * pi * cos (n*pi*x);
%! method_data.degree     = 3;     % Degree of the splines
%! method_data.regularity = 2;     % Regularity of the splines
%! method_data.nsub       = 20;    % Number of subdivisions
%! method_data.nquad      = 4;     % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u] = solve_laplace (problem_data, method_data);
%! [error_h1, error_l2] = sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex);
%! assert (msh.nel, 20)
%! assert (space.ndof, 23)
%! assert (error_h1, 0.824973349002851, 1e-14)
%! assert (error_l2, 0.00870630412762487, 1e-15)
