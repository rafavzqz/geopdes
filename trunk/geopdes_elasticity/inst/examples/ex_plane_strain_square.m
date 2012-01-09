% EX_PLANE_STRAIN_SQUARE: solve the plane-strain problem on a square.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = nrb4surf([0 0], [1 0], [0 1], [1 1]);

% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.press_sides  = [];
problem_data.drchlt_sides = [1 2 3 4];
problem_data.symm_sides   = [];

% Physical parameters
E  =  1; nu = .3; 
problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
fx = @(x, y) -(-(problem_data.mu_lame(x, y)*3 + problem_data.lambda_lame(x, y)).*sin(2*pi*x).*sin(2*pi*y) + ...
     (problem_data.mu_lame(x, y) + problem_data.lambda_lame(x, y)).*cos(2*pi*x).*cos(2*pi*y))*(2*pi)^2;
fy = fx;
problem_data.f = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));
problem_data.h       = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));

% Exact solution (optional)
uxex = @(x,y) sin(2*pi*x).*(sin(2*pi*y));
uyex = @(x,y) sin(2*pi*x).*(sin(2*pi*y));
problem_data.uex = @(x, y) cat(1, ...
                reshape (uxex (x,y), [1, size(x)]), ...
                reshape (uyex (x,y), [1, size(x)]));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3];     % Degree of the bsplines
method_data.regularity = [2 2];     % Regularity of the splines
method_data.nsub       = [9 9];     % Number of subdivisions
method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_plane_strain_2d (problem_data, method_data);

% 4) POST-PROCESSING. 
% 4.1) Export to Paraview
output_file = 'plane_strain_square_Deg3_Reg2_Sub9';

vtk_pts = {linspace(0, 1, 21), linspace(0, 1, 21)};
fprintf ('results being saved in: %s_displacement\n \n', output_file)
sp_to_vtk (u, space, geometry, vtk_pts, sprintf ('%s_displacement.vts', output_file), 'displacement')
sp_to_vtk_stress (u, space, geometry, vtk_pts, problem_data.lambda_lame, ...
                  problem_data.mu_lame, sprintf ('%s_stress', output_file)); 

% 4.2) Plot in Matlab. Comparison with the exact solution.
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

subplot (1,2,1)
quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
title ('Numerical solution'), axis equal tight
subplot (1,2,2)
eu2 = problem_data.uex (X, Y);
quiver (X, Y, squeeze(eu2(1,:,:)), squeeze(eu2(2,:,:)))
title ('Exact solution'), axis equal tight

error_l2 = sp_l2_error (space, msh, u, problem_data.uex)

%!demo
%! ex_plane_strain_square
