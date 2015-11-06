% EX_MAXWELL_SRC_CAMEMBERT: solve Maxwell source problem in three quarters of the cylinder, where the exact solution is a singular function.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_camembert.txt';

% Type of boundary conditions
problem_data.nmnn_sides   = [3 4 5 6];
problem_data.drchlt_sides = [1 2];

% Physical parameters
problem_data.c_mass  = @(x, y, z) ones(size(x));
problem_data.c_stiff = @(x, y, z) ones(size(x));

% Source and boundary terms
k = 1; % Constant that characterizes the singularity
problem_data.f = @(x, y, z) cat (1, ...
                          singular_function_maxwell (x, y, k), ...
                          zeros ([1, size(x)]));
problem_data.g = @(x, y, z, ind) zeros ([3, size(x)]);
problem_data.h = @(x, y, z, ind) cat (1, ...
                          singular_function_maxwell (x, y, k), ...
                          zeros ([1, size(x)]));

% Exact solution (optional)
problem_data.uex     = @(x, y, z) cat (1, ...
                          singular_function_maxwell (x, y, k), ...
                          zeros ([1, size(x)]));
problem_data.curluex = @(x, y, z) zeros ([3, size(x)]);

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [2 2 2];     % Degree of the bsplines
method_data.regularity = [1 1 1];     % Regularity of the splines
method_data.nsub       = [3 3 3];     % Number of subdivisions
method_data.nquad      = [3 3 3];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
% Remove warnings for the degenerated elements
warning ('off', 'geopdes:jacdet_zero_at_quad_node')
warning ('off', 'geopdes:zero_measure_element')
[geometry, msh, space, u] = solve_maxwell_src (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW
output_file = 'maxwell_camembert_Deg2_Reg1_Sub3';

vtk_pts = {linspace(0, 1, 15), linspace(0, 0.99, 15), linspace(0, 1, 15)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% 4.2) Comparison with the exact solution
[error_hcurl, error_l2] = ...
    sp_hcurl_error (space, msh, u, problem_data.uex, problem_data.curluex)

%!demo
%! ex_maxwell_src_camembert

%!test
%! problem_data.geo_name = 'geo_camembert.txt';
%! problem_data.nmnn_sides   = [3 4 5 6];
%! problem_data.drchlt_sides = [1 2];
%! problem_data.c_mass  = @(x, y, z) ones(size(x));
%! problem_data.c_stiff = @(x, y, z) ones(size(x));
%! k = 1; % Constant that characterizes the singularity
%! problem_data.f = @(x, y, z) cat (1, ...
%!                           singular_function_maxwell (x, y, k), ...
%!                           zeros ([1, size(x)]));
%! problem_data.g = @(x, y, z, ind) zeros ([3, size(x)]);
%! problem_data.h = @(x, y, z, ind) cat (1, ...
%!                           singular_function_maxwell (x, y, k), ...
%!                           zeros ([1, size(x)]));
%! problem_data.uex     = @(x, y, z) cat (1, ...
%!                           singular_function_maxwell (x, y, k), ...
%!                           zeros ([1, size(x)]));
%! problem_data.curluex = @(x, y, z) zeros ([3, size(x)]);
%! method_data.degree     = [2 2 2];     % Degree of the bsplines
%! method_data.regularity = [1 1 1];     % Regularity of the splines
%! method_data.nsub       = [3 3 3];     % Number of subdivisions
%! method_data.nquad      = [3 3 3];     % Points for the Gaussian quadrature rule
%! warning ('off', 'geopdes:jacdet_zero_at_quad_node')
%! warning ('off', 'geopdes:zero_measure_element')
%! [geometry, msh, space, u] = solve_maxwell_src (problem_data, method_data);
%! [error_hcurl, error_l2] = sp_hcurl_error (space, msh, u, problem_data.uex, problem_data.curluex);
%! assert (msh.nel, 81)
%! assert (space.ndof, 820)
%! assert (error_l2, 0.0436842163001183, 1e-14)
%! assert (error_hcurl, 0.0436861291289525, 1e-14)
