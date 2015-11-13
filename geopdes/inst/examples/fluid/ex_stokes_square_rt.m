% EX_STOKES_SQUARE_RT: solve the Stokes problem in the unit square with generalized Raviart-Thomas elements.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_square.txt';

% Type of boundary conditions for each side of the domain
problem_data.drchlt_sides = 1:4;
problem_data.nmnn_sides = [];

% Physical parameters
problem_data.viscosity = @(x, y) ones (size (x));

% Force term
fx = @(x, y) (6 * x + y .* cos(x .* y) + 2 * cos(y) .* sin(x));
fy = @(x, y) (x .* cos(x .* y) - 2 * cos(x) .* sin(y));
problem_data.f  = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));

% Boundary terms
uxex = @(x, y) (sin(x) .* cos(y));
uyex = @(x, y) (-sin(y) .* cos(x));

problem_data.h = @(x, y, iside) cat(1, ...
                  reshape (uxex (x,y), [1, size(x)]), ...
                  reshape (uyex (x,y), [1, size(x)]));

% Exact solution, to compute the errors
problem_data.velex = @(x, y) cat(1, ...
                  reshape (uxex (x,y), [1, size(x)]), ...
                  reshape (uyex (x,y), [1, size(x)]));

problem_data.pressex = @(x, y) (3 * x.^2 + sin(x .* y) - 1.239811742000564725943866);

problem_data.gradvelex = @test_stokes_square_bc_graduex;

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.element_name = 'rt';   % Element type for discretization
method_data.degree       = [ 3  3];  % Degree of the splines (pressure space)
method_data.regularity   = [ 2  2];  % Regularity of the splines (pressure space)
method_data.nsub         = [10 10];  % Number of subdivisions
method_data.nquad        = [ 5  5];  % Points for the Gaussian quadrature rule

% Penalization parameter for Nitsche's method
factor = 10;
method_data.Cpen = factor*(min(method_data.degree)+1);

% 3) CALL TO THE SOLVER
[geometry, msh, space_v, vel, space_p, press] = ...
                       solve_stokes (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) COMPARISON WITH EXACT SOLUTION
error_l2_p = sp_l2_error (space_p, msh, press, problem_data.pressex)
[error_h1_v, error_l2_v] = ...
   sp_h1_error (space_v, msh, vel, problem_data.velex, problem_data.gradvelex)


% 4.2) EXPORT TO PARAVIEW
output_file = 'SQUARE_RT_Deg3_Reg2_Sub10';

fprintf ('The result is saved in the files %s \n and %s \n \n', ...
           [output_file '_vel'], [output_file '_press']);
vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
sp_to_vtk (press, space_p, geometry, vtk_pts, [output_file '_press'], 'press')
sp_to_vtk (vel,   space_v, geometry, vtk_pts, [output_file '_vel'  ], {'velocity', 'divergence'}, {'value', 'divergence'})

% 4.3) PLOT IN MATLAB
[eu, F] = sp_eval (vel, space_v, geometry, vtk_pts);
[X,  Y] = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

figure()
subplot (1,2,2)
eu2 = problem_data.velex (X, Y);
quiver (X, Y, squeeze(eu2(1,:,:)), squeeze(eu2(2,:,:)))
axis equal
title('Exact solution')
subplot (1,2,1)
quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
axis equal
title('Computed solution')

%!demo
%! ex_stokes_square_rt

%!test
%! problem_data.geo_name = 'geo_square.txt';
%! problem_data.drchlt_sides = 1:4;
%! problem_data.nmnn_sides = [];
%! problem_data.viscosity = @(x, y) ones (size (x));
%! fx = @(x, y) (6 * x + y .* cos(x .* y) + 2 * cos(y) .* sin(x));
%! fy = @(x, y) (x .* cos(x .* y) - 2 * cos(x) .* sin(y));
%! problem_data.f  = @(x, y) cat(1, ...
%!                 reshape (fx (x,y), [1, size(x)]), ...
%!                 reshape (fy (x,y), [1, size(x)]));
%! uxex = @(x, y) (sin(x) .* cos(y));
%! uyex = @(x, y) (-sin(y) .* cos(x));
%! problem_data.h = @(x, y, iside) cat(1, ...
%!                   reshape (uxex (x,y), [1, size(x)]), ...
%!                   reshape (uyex (x,y), [1, size(x)]));
%! problem_data.velex = @(x, y) cat(1, ...
%!                   reshape (uxex (x,y), [1, size(x)]), ...
%!                   reshape (uyex (x,y), [1, size(x)]));
%! problem_data.pressex = @(x, y) (3 * x.^2 + sin(x .* y) - 1.239811742000564725943866);
%! problem_data.gradvelex = @test_stokes_square_bc_graduex;
%! method_data.element_name = 'rt';   % Element type for discretization
%! method_data.degree       = [ 3  3];  % Degree of the splines (pressure space)
%! method_data.regularity   = [ 2  2];  % Regularity of the splines (pressure space)
%! method_data.nsub         = [10 10];  % Number of subdivisions
%! method_data.nquad        = [ 5  5];  % Points for the Gaussian quadrature rule
%! factor = 10;
%! method_data.Cpen = factor*(min(method_data.degree)+1);
%! [geometry, msh, space_v, vel, space_p, press] = ...
%!                        solve_stokes (problem_data, method_data);
%! error_l2_p = sp_l2_error (space_p, msh, press, problem_data.pressex);
%! [error_h1_v, error_l2_v] = ...
%!    sp_h1_error (space_v, msh, vel, problem_data.velex, problem_data.gradvelex);
%! div = sp_eval (vel, space_v, geometry, [20 20], 'divergence');
%! assert (msh.nel, 100)
%! assert (space_p.ndof, 169)
%! assert (space_v.ndof, 364)
%! assert (error_l2_p, 4.93866840636771e-08, 1e-15)
%! assert (error_h1_v, 3.50559200757476e-06, 1e-15)
%! assert (error_l2_v, 5.53697801153355e-08, 1e-15)
%! assert (max (abs (div(:))) < 1e-12)