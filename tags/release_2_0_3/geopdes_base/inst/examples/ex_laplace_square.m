% EX_LAPLACE_SQUARE: solve the Poisson problem in the unit square with a B-spline discretization.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_square.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));

% Source and boundary terms
problem_data.f = @(x, y) zeros (size (x));
problem_data.g = @test_square_g_nmnn;
problem_data.h = @(x, y, ind) exp (x) .* sin(y);

% Exact solution (optional)
problem_data.uex     = @(x, y) exp (x) .* sin (y);
problem_data.graduex = @(x, y) cat (1, ...
                       reshape (exp(x).*sin(y), [1, size(x)]), ...
                       reshape (exp(x).*cos(y), [1, size(x)]));

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

output_file = 'Square_BSP_Deg3_Reg2_Sub9';

vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% 4.2) PLOT IN MATLAB. COMPARISON WITH THE EXACT SOLUTION

[eu, F] = sp_eval (u, space, geometry, vtk_pts);
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
%! ex_laplace_square
