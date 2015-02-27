% EX_LAPLACE_CUBE: solve the Poisson problem in the unit cube with a B-spline discretization.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_cube.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [4 5 6];
problem_data.drchlt_sides = [1 2 3];

% Physical parameters
problem_data.c_diff  = @(x, y, z) ones(size(x));

% Source and boundary terms
problem_data.f = @(x, y, z) -exp(x+z).*sin(y);
problem_data.g = @test_cube_g_nmnn;
problem_data.h = @(x, y, z, ind) exp (x + z) .* sin (y);

% Exact solution (optional)
problem_data.uex     = @(x, y, z) exp (x + z) .* sin (y);
problem_data.graduex = @(x, y, z) cat (1, ...
                          reshape (exp (x + z) .* sin (y), [1, size(x)]), ...
                          reshape (exp (x + z) .* cos (y), [1, size(x)]), ...
                          reshape (exp (x + z) .* sin (y), [1, size(x)]));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3 3];       % Degree of the splines
method_data.regularity = [2 2 2];       % Regularity of the splines
method_data.nsub       = [3 3 3];       % Number of subdivisions
method_data.nquad      = [4 4 4];       % Points for Gaussian quadrature rule

% 3) CALL TO THE SOLVER

[geometry, msh, space, u] = solve_laplace_3d (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW

output_file = 'Cube_BSP_Deg3_Reg2_Sub3';

vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20), linspace(0, 1, 20)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% 4.2) COMPARISON WITH THE EXACT SOLUTION

[error_h1, error_l2] = ...
           sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex)

%!demo
%! ex_laplace_cube
