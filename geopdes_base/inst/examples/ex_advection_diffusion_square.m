% EX_ADVECTION_DIFFUSION_SQUARE: solve the advection-diffusion problem in the unit square.

% PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_square.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];

% Physical parameters
problem_data.c_diff  = @(x, y) 1e-3 * ones(size(x)); % diffusion coefficient (mu)
problem_data.grad_diff = @(x, y) zeros ([2, size(x)]);  % grad(mu)
theta = pi / 4;
problem_data.vel = @(x, y) cat (1, ...
             reshape (-cos(theta) * ones(size(x)), [1, size(x)]), ...
             reshape (-sin(theta) * ones(size(x)), [1, size(x)]));

problem_data.f = @(x, y) zeros (size (x));
problem_data.g = @(x, y, ind) zeros (size (x));
problem_data.h = @test_adv_diff_h_drchlt;

% CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3];    % Degree of the splines
method_data.regularity = [2 2];    % Regularity of the splines
method_data.nsub       = [20 20];  % Number of subdivisions
method_data.nquad      = [4 4];    % Points for the Gaussian quadrature rule
method_data.stab = true; % Choose whether to perform stabilization (true/false)

% CALL TO THE SOLVER
[geometry, msh, space, u] = solve_adv_diff_2d (problem_data, method_data);

% POST-PROCESSING
% Export to Paraview

output_file = 'Advection_Diffusion_square';

vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% Plot in Matlab

[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
surf (X, Y, eu)

disp ('The overshoots come from using the L2 projection for Dirichlet boundary conditions')
disp ('This can be fixed setting the boundary conditions in a different way (see solve_adv_diff_2d)')


%!demo
%! ex_advection_diffusion_square
