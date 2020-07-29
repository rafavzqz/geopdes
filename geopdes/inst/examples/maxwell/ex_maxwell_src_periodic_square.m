% EX_MAXWELL_SRC_PERIODIC_SQUARE: solve Maxwell source problem in the unit square.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_square.txt';

% Type of boundary conditions
problem_data.nmnn_sides     = [];
problem_data.drchlt_sides   = [3 4];
problem_data.periodic_sides = [1 2];

% Physical parameters
problem_data.c_stiff = @(x, y) ones(size(x));
problem_data.c_mass  = @(x, y) pi^2 * ones(size(x));

% Source and boundary terms
problem_data.f = @(x, y) cat(1, zeros (1, size (x, 1), size (x, 2)), ...
                             reshape (5*pi^2*sin(2*pi*x), [1, size(x)]));
problem_data.g = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
problem_data.h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));

% Exact solution (optional)
problem_data.uex     = @(x, y) cat(1, zeros (1, size (x, 1), size (x, 2)), ...
                                      reshape ( sin(2*pi*x), [1, size(x)]));
problem_data.curluex = @(x, y) 2*pi*cos(2*pi*x);

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data 
method_data.degree     = [3 3];     % Degree of the bsplines
method_data.regularity = [2 2];     % Regularity of the splines
method_data.nsub       = [16 16];     % Number of subdivisions
method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_maxwell_src (problem_data, method_data);

% % % % 4) POST-PROCESSING
% % % vtk_pts = {linspace(0, 1, 30), linspace(0, 1, 30)};
% % % % 4.1) EXPORT TO PARAVIEW
% % % output_file = 'maxwell_square_Deg3_Reg2_Sub8';
% % % fprintf ('The result is saved in the file %s \n \n', output_file);
% % % sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% 4.2) Plot in Matlab. Comparison with the exact solution
% % % [eu, F] = sp_eval (u, space, geometry, vtk_pts);
% % % [X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
% % % eu2     = problem_data.uex (X, Y);
% % % 
% % % subplot(1,2,1)
% % % quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
% % % axis equal tight
% % % title('Computed solution')
% % % ylim([0,1]);
% % % subplot(1,2,2)
% % % quiver (X, Y, squeeze(eu2(1,:,:)), squeeze(eu2(2,:,:)))
% % % axis equal tight
% % % title('Exact solution')
% % % ylim([0,1]);


y_pt = 0.75;
vtk_pts = {linspace(0,1,30),y_pt};
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
eu2     = problem_data.uex (X, Y);

eu = squeeze(eu(2,:,:)); %eu = eu(:,1);
eu2 = squeeze(eu2(2,:,:)); %eu2 = eu2(:,1);

subplot(1,2,1)
plot (X, eu);
grid on;
title('Computed solution')

subplot(1,2,2)
plot (X, eu2)
grid on;
title('Exact solution')

[error_hcurl, error_l2] = ...
    sp_hcurl_error (space, msh, u, problem_data.uex, problem_data.curluex)

%!demo
%! ex_maxwell_src_square

%!test
%! problem_data.geo_name = 'geo_square.txt';
%! problem_data.nmnn_sides   = [2 3];
%! problem_data.drchlt_sides = [1 4];
%! problem_data.c_stiff = @(x, y) ones(size(x));
%! problem_data.c_mass  = @(x, y) ones(size(x));
%! problem_data.f = @(x, y) cat(1, ...
%!                reshape (-exp(x) .* sin(y) + 2*sin(y), [1, size(x)]), ...
%!                zeros ([1, size(x)]));
%! problem_data.g = @(x, y, ind) test_maxwell_square_g_nmnn (x, y, ind);
%! problem_data.h = @(x, y, ind) cat(1, ...
%!                     reshape (sin(y), [1, size(x)]), ...
%!                     reshape (exp(x) .* cos(y), [1, size(x)]));
%! problem_data.uex     = @(x, y) cat(1, ...
%!                     reshape (sin(y), [1, size(x)]), ...
%!                     reshape (exp(x) .* cos(y), [1, size(x)]));
%! problem_data.curluex = @(x, y) exp(x).*cos(y) - cos(y);
%! method_data.degree     = [3 3];     % Degree of the bsplines
%! method_data.regularity = [2 2];     % Regularity of the splines
%! method_data.nsub       = [8 8];     % Number of subdivisions
%! method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u] = solve_maxwell_src (problem_data, method_data);
%! [error_hcurl, error_l2] = sp_hcurl_error (space, msh, u, problem_data.uex, problem_data.curluex);
%! assert (msh.nel, 64)
%! assert (space.ndof, 220)
%! assert (error_l2, 9.90499820528337e-06, 1e-16)
%! assert (error_hcurl, 1.94621849955079e-05, 1e-16)