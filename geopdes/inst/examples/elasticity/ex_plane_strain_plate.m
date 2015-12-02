% EX_PLANE_STRAIN_PLATE: solve the plane-strain problem on a square plate with a circular hole.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_plate_with_hole.txt';

% Type of boundary conditions
problem_data.nmnn_sides   = [4];
problem_data.drchlt_sides = [];
problem_data.press_sides  = [];
problem_data.symm_sides   = [1 2];

% Physical parameters
E  =  1e5; nu = .3; 
problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
problem_data.f = @(x, y) zeros (2, size (x, 1), size (x, 2));
problem_data.g = @test_plate_with_hole_g_nmnn;
problem_data.h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3];     % Degree of the basis functions
method_data.regularity = [2 2];     % Regularity of the basis functions
method_data.nsub       = [8 8];     % Number of subdivisions
method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_linear_elasticity (problem_data, method_data);

% 4) POST-PROCESSING. 
% 4.1) Export to Paraview
output_file = 'plane_strain_plate_Deg3_Reg2_Sub8';

vtk_pts = {linspace(0, 1, 21), linspace(0, 1, 21)};
fprintf ('results being saved in: %s \n \n', output_file)
sp_to_vtk (u, space, geometry, vtk_pts, output_file, {'displacement', 'stress'}, {'value', 'stress'}, ...
    problem_data.lambda_lame, problem_data.mu_lame)

% 4.2) Plot in Matlab
[eu, F] = sp_eval (u, space, geometry, vtk_pts, {'value', 'stress'}, problem_data.lambda_lame, problem_data.mu_lame);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

figure
quiver (X, Y, squeeze(eu{1}(1,:,:)), squeeze(eu{1}(2,:,:)))
axis equal tight
title ('Computed displacement')

figure
surf (X, Y, squeeze(eu{2}(1,1,:,:)))
view(2); axis equal; shading interp; colorbar;
title ('Stress component \sigma_{xx}')

%!demo
%! ex_plane_strain_plate
