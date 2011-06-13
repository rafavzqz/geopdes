% EX_LAPLACE_RING_MIXED_BC: solve the Poisson problem in one quarter of a ring, discretized with B-splines (non-isoparametric approach).

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_ring.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [1 3 4];
problem_data.drchlt_sides = [2];

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));

% Source and boundary terms
problem_data.f = @(x, y) exp(x).*((x.^2 + y.^2 - 1).*sin(x.*y) - 2*y.*cos(x.*y));
problem_data.g = @test_ring_mixed_bc_g_nmnn;
problem_data.h = @(x, y, ind) exp(x).*sin(x.*y);

% Exact solution (optional)
problem_data.uex     = @(x, y) exp(x) .* sin (x.*y);
problem_data.graduex = @(x, y) cat (1, ...
               reshape ( exp(x).*(sin(x.*y) + y.*cos(x.*y)), [1, size(x)]), ...
               reshape (exp(x).*x.*cos(x.*y), [1, size(x)]));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data 
method_data.degree     = [3 3];       % Degree of the splines
method_data.regularity = [2 2];       % Regularity of the splines
method_data.nsub       = [9 9];       % Number of subdivisions
method_data.nquad      = [4 4];       % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_laplace_2d (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW

output_file = 'Ring_mixed_bc_BSP_Deg3_Reg2_Sub9';

vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk_new (u, space, geometry, vtk_pts, output_file, 'u')

% 4.2) PLOT IN MATLAB. COMPARISON WITH THE EXACT SOLUTION

[eu, F] = sp_eval_pts (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
subplot (1,2,1)
surf (X, Y, eu)
title ('Numerical solution'), axis tight
subplot (1,2,2)
surf (X, Y, problem_data.uex (X,Y))
title ('Exact solution'), axis tight

% Display errors of the computed solution in the L2 and H1 norm
[error_h1, error_l2] = ...
           sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex)

%!demo
%! ex_laplace_ring_mixed_bc
