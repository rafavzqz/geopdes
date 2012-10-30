% EX_LAPLACE_RING_TSPLINES: solve the Poisson problem in one quarter of a ring
%  with a T-spline discretization, coming from the T-splines plug-in for Rhino.

% PHYSICAL DATA OF THE PROBLEM
clear problem_data
% The physical domain is one quarter of a ring, defined either as a NURBS (ring_example) 
%  or as a structured T-spline with T-junctions (arc_tsplines).
% problem_data.geo_name = 'arch_nurbs.iga';
problem_data.geo_name = 'arch_tsplines.iga';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [1 3 4];
problem_data.drchlt_sides = [2];

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));

% Source and boundary terms
problem_data.f = @(x, y) exp(x).*((x.^2 + y.^2 - 1).*sin(x.*y) - 2*y.*cos(x.*y));
problem_data.g = @test_ring_mixed_bc_g_nmnn;
problem_data.h = @(x, y, ind) exp(x) .* sin(x.*y);

% Exact solution (optional)
problem_data.uex     = @(x, y) exp(x) .* sin (x.*y);
problem_data.graduex = @(x, y) cat (1, ...
               reshape (exp(x).*(sin(x.*y) + y.*cos(x.*y)), [1, size(x)]), ...
               reshape (exp(x).*x.*cos(x.*y), [1, size(x)]));

% CALL TO THE SOLVER
[tspline, msh, space, u] = solve_laplace_tsplines_2d (problem_data);

% POST-PROCESSING
% Display errors of the computed solution in the L2 and H1 norm
[error_h1, error_l2] = ...
           tspline_h1_error (space, msh, u, problem_data.uex, problem_data.graduex)

% Export to Paraview
output_file = 'Ring_tsplines';
vtk_pts_per_element = [10 10]; % Number of points on each Bezier element
tspline_to_vtk (u, tspline, vtk_pts_per_element, output_file, 'u');
