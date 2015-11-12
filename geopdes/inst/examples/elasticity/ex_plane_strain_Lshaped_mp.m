% EX_PLANE_STRAIN_LSHAPED_MP: solve plane strain problem in the multipatch L-shaped domain.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_Lshaped_mp.txt';
% problem_data.geo_name = 'geo_Lshaped_mp_b.txt';
% problem_data.geo_name = 'geo_Lshaped_mp_c.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4 5 6];

% Physical parameters
E  =  1; nu = 0.3; 
problem_data.lambda_lame = @(x, y, z) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x))); 
problem_data.mu_lame = @(x, y, z) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
fx = @(x, y) problem_data.lambda_lame (x, y) .* sin(y) + ...
             problem_data.mu_lame(x, y) .* (1 + x) .* sin(y);
fy = @(x, y) problem_data.mu_lame (x, y) .* cos(y) .* (2*x - 1) - ...
             problem_data.lambda_lame (x, y) .* cos(y) .* (1 - x);
problem_data.f  = @(x, y) cat(1, ...
                   reshape (fx (x,y), [1, size(x)]), ...
                   reshape (fy (x,y), [1, size(x)]));
problem_data.h = @(x, y, ind) cat(1, ...
                   reshape (x.*sin(y), [1, size(x)]), ...
                   reshape (x.*cos(y), [1, size(x)]));

% Exact solution (optional)
uxex = @(x, y) x.*sin(y);
uyex = @(x, y) x.*cos(y);
problem_data.uex  = @(x, y) cat(1, ...
		      reshape (uxex (x,y), [1, size(x)]), ...
		      reshape (uyex (x,y), [1, size(x)]));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3];     % Degree of the bsplines
method_data.regularity = [2 2];     % Regularity of the splines
method_data.nsub       = [6 6];     % Number of subdivisions
method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = mp_solve_linear_elasticity (problem_data, method_data);

% 4) POST-PROCESSING
% EXPORT TO PARAVIEW
output_file = 'plane_strain_Lshaped_mp_Deg3_Reg2_Sub6';

vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
fprintf ('results being saved in: %s.pvd\n \n', output_file)
sp_to_vtk (u, space, geometry, vtk_pts, output_file, {'displacement', 'stress'}, {'value', 'stress'}, ...
    problem_data.lambda_lame, problem_data.mu_lame)

% COMPARISON WITH THE EXACT SOLUTION
error_l2 = sp_l2_error (space, msh, u, problem_data.uex)

%!demo
%! ex_plane_strain_Lshaped_mp

%!test
%! problem_data.geo_name = 'geo_Lshaped_mp.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! E  =  1; nu = 0.3; 
%! problem_data.lambda_lame = @(x, y, z) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x))); 
%! problem_data.mu_lame = @(x, y, z) (E/(2*(1+nu)) * ones (size (x)));
%! fx = @(x, y) problem_data.lambda_lame (x, y) .* sin(y) + ...
%!              problem_data.mu_lame(x, y) .* (1 + x) .* sin(y);
%! fy = @(x, y) problem_data.mu_lame (x, y) .* cos(y) .* (2*x - 1) - ...
%!              problem_data.lambda_lame (x, y) .* cos(y) .* (1 - x);
%! problem_data.f  = @(x, y) cat(1, ...
%!                    reshape (fx (x,y), [1, size(x)]), ...
%!                    reshape (fy (x,y), [1, size(x)]));
%! problem_data.h = @(x, y, ind) cat(1, ...
%!                         reshape (x.*sin(y), [1, size(x)]), ...
%! 			reshape (x.*cos(y), [1, size(x)]));
%! uxex = @(x, y) x.*sin(y);
%! uyex = @(x, y) x.*cos(y);
%! problem_data.uex  = @(x, y) cat(1, ...
%! 		      reshape (uxex (x,y), [1, size(x)]), ...
%! 		      reshape (uyex (x,y), [1, size(x)]));
%! method_data.degree     = [3 3];     % Degree of the bsplines
%! method_data.regularity = [2 2];     % Regularity of the splines
%! method_data.nsub       = [6 6];     % Number of subdivisions
%! method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u] = mp_solve_linear_elasticity (problem_data, method_data);
%! error_l2 = sp_l2_error (space, msh, u, problem_data.uex);
%! assert (error_l2, 6.41633220768243e-07, 1e-16)
%! assert (space.ndof, 450);
%!
%! problem_data.geo_name = 'geo_Lshaped_mp_b.txt';
%! [geometry, msh, space, u] = mp_solve_linear_elasticity (problem_data, method_data);
%! error_l2 = sp_l2_error (space, msh, u, problem_data.uex);
%! assert (error_l2, 6.41633220768243e-07, 1e-16)
%! assert (space.ndof, 450);
%!
%! problem_data.geo_name = 'geo_Lshaped_mp_c.txt';
%! [geometry, msh, space, u] = mp_solve_linear_elasticity (problem_data, method_data);
%! error_l2 = sp_l2_error (space, msh, u, problem_data.uex);
%! assert (error_l2, 6.41633220768243e-07, 1e-16)
%! assert (space.ndof, 450);
