% EX_MAXWELL_SRC_LSHAPED: solve Maxwell source problem in the Lshaped domain, where the exact solution is a singular function.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_Lshaped_C0.txt';
%problem_data.geo_name = 'geo_Lshaped_C1.txt';

% Type of boundary conditions
problem_data.nmnn_sides   = [1 3 4];
problem_data.drchlt_sides = [2];

% Physical parameters
problem_data.c_mass  = @(x, y) ones(size(x));
problem_data.c_stiff = @(x, y) ones(size(x));

% Source and boundary terms
k = 1; % Constant that characterizes the singularity
problem_data.f = @(x, y) singular_function_maxwell (x, y, k);
problem_data.g = @(x, y, ind) zeros([2, size(x)]);
problem_data.h = @(x, y, ind) singular_function_maxwell (x, y, k);

% Exact solution (optional)
problem_data.uex     = @(x, y) singular_function_maxwell (x, y, k);
problem_data.curluex = @(x, y) zeros (size(x));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data 
method_data.degree     = [3 3];     % Degree of the bsplines
method_data.regularity = [2 2];     % Regularity of the splines
method_data.nsub       = [8 8];     % Number of subdivisions
method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_maxwell_src (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW
output_file = 'maxwell_Lshaped_Deg3_Reg2_Sub8';

vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% 4.2) Plot in Matlab. Comparison with the exact solution
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
eu2     = problem_data.uex (X, Y);

subplot(1,2,1)
quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
axis equal tight
title('Computed solution')
subplot(1,2,2)
quiver (X, Y, squeeze(eu2(1,:,:)), squeeze(eu2(2,:,:)))
axis equal tight
title('Exact solution')

[error_hcurl, error_l2] = ...
    sp_hcurl_error (space, msh, u, problem_data.uex, problem_data.curluex)

%!demo
%! ex_maxwell_src_Lshaped

%!test
%! problem_data.geo_name = 'geo_Lshaped_C0.txt';
%! problem_data.nmnn_sides   = [1 3 4];
%! problem_data.drchlt_sides = [2];
%! problem_data.c_mass  = @(x, y) ones(size(x));
%! problem_data.c_stiff = @(x, y) ones(size(x));
%! k = 1; % Constant that characterizes the singularity
%! problem_data.f = @(x, y) singular_function_maxwell (x, y, k);
%! problem_data.g = @(x, y, ind) zeros([2, size(x)]);
%! problem_data.h = @(x, y, ind) singular_function_maxwell (x, y, k);
%! problem_data.uex     = @(x, y) singular_function_maxwell (x, y, k);
%! problem_data.curluex = @(x, y) zeros (size(x));
%! method_data.degree     = [3 3];     % Degree of the bsplines
%! method_data.regularity = [2 2];     % Regularity of the splines
%! method_data.nsub       = [8 8];     % Number of subdivisions
%! method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u] = solve_maxwell_src (problem_data, method_data);
%! [error_hcurl, error_l2] = sp_hcurl_error (space, msh, u, problem_data.uex, problem_data.curluex);
%! assert (msh.nel, 128)
%! assert (space.ndof, 430)
%! assert (error_l2, 0.0467860957084692, 1e-15)
%! assert (error_hcurl, 0.0467926161512472, 1e-15)
