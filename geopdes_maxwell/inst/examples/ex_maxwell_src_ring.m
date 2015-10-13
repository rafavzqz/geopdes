% EX_MAXWELL_SRC_RING: solve Maxwell source problem in one quarter of a ring.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_ring.txt';

% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];

% Physical parameters
problem_data.c_stiff = @(x, y) ones(size(x));
problem_data.c_mass  = @(x, y) ones(size(x));

% Source and boundary terms
problem_data.f = @(x, y) cat(1, ...
             reshape (2*sin(y), [1, size(x)]), ...
             reshape (2*sin(x), [1, size(x)]));
problem_data.g = @(x, y, ind) test_maxwell_ring_g_nmnn (x, y, ind);
problem_data.h = @(x, y, ind) cat(1, ...
                    reshape (sin(y), [1, size(x)]), ...
                    reshape (sin(x), [1, size(x)]));

% Exact solution (optional)
problem_data.uex     = @(x, y) cat(1, ...
                    reshape (sin(y), [1, size(x)]), ...
                    reshape (sin(x), [1, size(x)]));
problem_data.curluex = @(x, y) cos(x) - cos(y);

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data 
method_data.degree     = [3 3];     % Degree of the bsplines
method_data.regularity = [2 2];     % Regularity of the splines
method_data.nsub       = [9 9];     % Number of subdivisions
method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_maxwell_src (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW
output_file = 'maxwell_ring_Deg3_Reg2_Sub9';

vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% 4.2) Plot in Matlab. Comparison with the exact solution
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
eu2     = problem_data.uex (X, Y);

subplot(1,2,1),
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
%! ex_maxwell_src_ring

%!test
%! problem_data.geo_name = 'geo_ring.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4];
%! problem_data.c_stiff = @(x, y) ones(size(x));
%! problem_data.c_mass  = @(x, y) ones(size(x));
%! problem_data.f = @(x, y) cat(1, ...
%!              reshape (2*sin(y), [1, size(x)]), ...
%!              reshape (2*sin(x), [1, size(x)]));
%! problem_data.g = @(x, y, ind) test_maxwell_ring_g_nmnn (x, y, ind);
%! problem_data.h = @(x, y, ind) cat(1, ...
%!                     reshape (sin(y), [1, size(x)]), ...
%!                     reshape (sin(x), [1, size(x)]));
%! problem_data.uex     = @(x, y) cat(1, ...
%!                     reshape (sin(y), [1, size(x)]), ...
%!                     reshape (sin(x), [1, size(x)]));
%! problem_data.curluex = @(x, y) cos(x) - cos(y);
%! method_data.degree     = [3 3];     % Degree of the bsplines
%! method_data.regularity = [2 2];     % Regularity of the splines
%! method_data.nsub       = [9 9];     % Number of subdivisions
%! method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u] = solve_maxwell_src (problem_data, method_data);
%! [error_hcurl, error_l2] = sp_hcurl_error (space, msh, u, problem_data.uex, problem_data.curluex);
%! assert (msh.nel, 81)
%! assert (space.ndof, 264)
%! assert (error_l2, 3.85922857182716e-04, 1e-16)
%! assert (error_hcurl, 5.72616592482421e-04, 1e-16)
