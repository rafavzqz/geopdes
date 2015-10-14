% EX_KIRCHHOFF_CIRCULAR_PLATE: solve bilaplacian in a simply supported circular plate.

% PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Geometry definition (unrefined geometry)
Radius = 1.0;     % Radius of the plate

% Knot vector, control points (in homogeneous coordinates) and weights
knots = {[0 0 0 1 1 1], [0 0 0 1 1 1]};
coefs = zeros (4, 3, 3);
coefs(:, 1, 1) = [1, 0, 0, 1];
coefs(:, 2, 1) = [1/sqrt(2), -1/sqrt(2), 0 , 1/sqrt(2)];
coefs(:, 3, 1) = [0, -1, 0, 1];
coefs(:, 1, 2) = [1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
coefs(:, 2, 2) = [0, 0, 0, sqrt(2)-1];
coefs(:, 3, 2) = [-1/sqrt(2), -1/sqrt(2), 0, 1/sqrt(2)];
coefs(:, 1, 3) = [0, 1, 0, 1];
coefs(:, 2, 3) = [-1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
coefs(:, 3, 3) = [-1, 0, 0, 1];
srf = nrbmak (coefs, knots);
srf = nrbtform (srf, vecscale ([Radius Radius 1]));

problem_data.geo_name = srf;

% BOUNDARY CONDITIONS 
problem_data.simply_supported_sides = [1 2 3 4];
problem_data.clamped_sides = [];

% Physical parameters
D = 1;  % Flexural rigidity of the plate
problem_data.c_diff  = @(x, y) D*ones(size(x));

% Source term
p = -1;  % Distributed load
problem_data.f = @(x, y) p*ones(size(x));

% CHOICE OF THE DISCRETIZATION PARAMETERS
%clear method_data
method_data.degree     = [3 3];  % Degree of the splines
method_data.regularity = [2 2];  % Regularity of the splines
method_data.nsub       = [9 9];  % Number of subdivisions
method_data.nquad      = [4 4];  % Points for the Gaussian quadrature rule

% CALL TO THE SOLVER
[geometry, msh, space, u] =...
    solve_bilaplace_gradgrad_2d_iso (problem_data, method_data);

% POST-PROCESSING
% EXPORT TO PARAVIEW
output_file = 'Kirchhoff_Circular_Plate_Simply_Supported';
vtk_pts = {linspace(0, 1, 41), linspace(0, 1, 41)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% % PLOT IN MATLAB
% % Plot of the refined geometry and control points, and of the computational mesh
% figure
% nrbctrlplot (geometry.nurbs)
% figure
% nrbkntplot (geometry.nurbs)

% Plot of the computed solution
figure
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
surf (X, Y, eu)
title ('Numerical solution')
axis equal

% Max Deflection (in the same points used for plotting)
max_displacement = min (eu(:));
exact_solution = -5*Radius^4/(64*D);
fprintf('Computed solution, max. displacement = %e \n', max_displacement);
fprintf('Exact solution, max. displacement = %e \n', exact_solution);
