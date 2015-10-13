% EX_LIN_ELAST_HORSESHOE: solve the linear elasticity problem in the horseshoe.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'horseshoe.mat';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [1 2 3 4];
problem_data.drchlt_sides = [5 6];
problem_data.press_sides  = [];
problem_data.symm_sides   = [];

% Physical parameters
E  =  1; nu = 0.3;
problem_data.lambda_lame = @(x, y, z) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x))); 
problem_data.mu_lame = @(x, y, z) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
fx = @(x, y, z) zeros (size (x));
fy = fx;
fz = @(x, y, z) 0.1 * ones (size (x));
problem_data.f = @(x, y, z) cat(1, ...
                   reshape (fx (x,y,z), [1, size(x)]), ...
                   reshape (fy (x,y,z), [1, size(x)]), ...
                   reshape (fz (x,y,z), [1, size(x)]));
problem_data.g = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2));
problem_data.h = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3 3];     % Degree of the bsplines
method_data.regularity = [2 2 2];     % Regularity of the splines
method_data.nsub       = [1 1 1];     % Number of subdivisions
method_data.nquad      = [4 4 4];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_linear_elasticity (problem_data, method_data);

% 4) POST-PROCESSING. EXPORT TO PARAVIEW
output_file = 'lin_elast_horseshoe_Deg3_Reg2_Sub1';

vtk_pts = {linspace(0, 1, 8), linspace(0, 1, 8), linspace(0, 1, 40)};
fprintf ('results being saved in: %s \n \n', output_file)
sp_to_vtk (u, space, geometry, vtk_pts, output_file, {'displacement', 'stress'}, {'value', 'stress'}, ...
    problem_data.lambda_lame, problem_data.mu_lame)

% Plot in Matlab
figure
def_geom = geo_deform (u, space, geometry);
nrbplot (def_geom.nurbs, [40 40 40], 'light', 'on')
axis equal tight
title ('Deformed configuration')

%!demo
%! ex_lin_elast_horseshoe

%!test
%! problem_data.geo_name = 'horseshoe.mat';
%! problem_data.nmnn_sides   = [1 2 3 4];
%! problem_data.drchlt_sides = [5 6];
%! problem_data.press_sides  = [];
%! problem_data.symm_sides   = [];
%! E  =  1; nu = 0.3;
%! problem_data.lambda_lame = @(x, y, z) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x))); 
%! problem_data.mu_lame = @(x, y, z) (E/(2*(1+nu)) * ones (size (x)));
%! fx = @(x, y, z) zeros (size (x));
%! fy = fx;
%! fz = @(x, y, z) 0.1 * ones (size (x));
%! problem_data.f = @(x, y, z) cat(1, ...
%!                    reshape (fx (x,y,z), [1, size(x)]), ...
%!                    reshape (fy (x,y,z), [1, size(x)]), ...
%!                    reshape (fz (x,y,z), [1, size(x)]));
%! problem_data.g = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2));
%! problem_data.h = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2));
%! method_data.degree     = [3 3 3];     % Degree of the bsplines
%! method_data.regularity = [2 2 2];     % Regularity of the splines
%! method_data.nsub       = [1 1 1];     % Number of subdivisions
%! method_data.nquad      = [4 4 4];     % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u] = solve_linear_elasticity (problem_data, method_data);
%! assert (msh.nel, 12)
%! assert (space.ndof, 1080)