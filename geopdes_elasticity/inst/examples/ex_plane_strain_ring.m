% EX_PLANE_STRAIN_RING: solve the plane-strain problem on one quarter of a cylinder.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_ring.txt';

% Type of boundary conditions
problem_data.nmnn_sides   = [2];
problem_data.drchlt_sides = [];
problem_data.press_sides  = [1];
problem_data.symm_sides   = [3 4];

% Physical parameters
E  =  1; nu = 0; 
problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
P = 1;
problem_data.f = @(x, y) zeros (2, size (x, 1), size (x, 2));
problem_data.g = @(x, y, ind) test_plane_strain_ring_g_nmnn (x, y, P, nu, ind);
problem_data.h = @(x, y, ind) test_plane_strain_ring_uex (x, y, E, nu, P);
problem_data.p = @(x, y, ind) P * ones (size (x));

% Exact solution (optional)
problem_data.uex = @(x, y) test_plane_strain_ring_uex (x, y, E, nu, P);

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3];     % Degree of the basis functions
method_data.regularity = [2 2];     % Regularity of the basis functions
method_data.nsub       = [9 9];     % Number of subdivisions
method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_plane_strain_2d (problem_data, method_data);

% 4) POST-PROCESSING. 
% 4.1) Export to Paraview
output_file = 'plane_strain_ring_Deg3_Reg2_Sub9';

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
%! ex_plane_strain_ring
