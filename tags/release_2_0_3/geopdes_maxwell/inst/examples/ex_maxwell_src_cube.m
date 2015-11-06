% EX_MAXWELL_SRC_CUBE: solve Maxwell source problem in the unit cube.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_cube.txt';

% Type of boundary conditions
problem_data.nmnn_sides   = [3 4 5 6];
problem_data.drchlt_sides = [1 2];

% Physical parameters
problem_data.c_mass  = @(x, y, z) ones(size(x));
problem_data.c_stiff = @(x, y, z) ones(size(x));

% Source and boundary terms
problem_data.f = @(x, y, z) cat(1, ...
                    reshape (sin(y) .* (2*z - exp(x)) - exp(z) .* sin(x), [1, size(x)]), ...
                    zeros ([1, size(x)]), ...
                    reshape (2 * exp(z) .* cos(x), [1, size(x)]));
problem_data.g = @(x, y, z, ind) test_maxwell_cube_g_nmnn (x, y, z, ind);
problem_data.h = @(x, y, z, ind) test_maxwell_cube_h_drchlt (x, y, z, ind);

% Exact solution (optional)
problem_data.uex     = @(x, y, z) cat(1, ...
                    reshape (sin(y) .* z, [1, size(x)]), ...
                    reshape (exp(x) .* cos(y), [1, size(x)]), ...
                    reshape (exp(z) .* cos(x), [1, size(x)]));
problem_data.curluex = @(x, y, z) cat(1, ...
                    zeros ([1, size(x)]), ...
                    reshape (sin(y) + exp(z) .* sin(x), [1, size(x)]), ...
                    reshape (exp(x).*cos(y) - z .* cos(y), [1, size(x)]));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [2 2 2];     % Degree of the bsplines
method_data.regularity = [1 1 1];     % Regularity of the splines
method_data.nsub       = [3 3 3];     % Number of subdivisions
method_data.nquad      = [3 3 3];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_maxwell_src_3d (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW
output_file = 'maxwell_cube_Deg2_Reg1_Sub3';

vtk_pts = {linspace(0, 1, 15), linspace(0, 1, 15), linspace(0, 1, 15)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% 4.2) Comparison with the exact solution
[error_hcurl, error_l2] = ...
    sp_hcurl_error (space, msh, u, problem_data.uex, problem_data.curluex)

%!demo
%! ex_maxwell_src_cube
