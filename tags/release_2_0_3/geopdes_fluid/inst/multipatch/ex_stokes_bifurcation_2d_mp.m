% EX_STOKES_BIFURCATION_MP: data file for Stokes problem in a pipe with a bifurcation.

% 1) PHYSICAL DATA OF THE PROBLEM
problem_data  = struct ();

% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_bifurcation_mp.txt';


% Type of boundary conditions for each side of the domain
problem_data.drchlt_sides = 1:3;
problem_data.nmnn_sides = [];

% Physical parameters
problem_data.viscosity = @(x, y, z) ones (size (x));

% Force term
problem_data.f  = @(x, y, z) zeros ([2, size(x)]);

% Boundary terms
problem_data.h  = @test_stokes_bifurcation_mp_h_drchlt;

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
method_data = struct ();

method_data.element_name = 'th';     % Element type for discretization

method_data.degree       = [3 3];  % Degree of the splines
method_data.regularity   = [2 2];  % Regularity of the splines
method_data.nsub         = [5 5];  % Number of subdivisions
method_data.nquad        = [5 5];  % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space_v, vel, gnum, space_p, press, gnump] = mp_solve_stokes_2d (problem_data, method_data);


% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW
output_file  = 'bifurcation_2d_mp_deg3_reg2_sub5';
vtk_pts = {linspace(0, 1, 10), linspace(0, 1, 10)};

fprintf ('results being saved in: %s_vel.pvd and %s_press.pvd\n', output_file, output_file)
mp_sp_to_vtk (vel, space_v, geometry, gnum, vtk_pts, sprintf ('%s_vel', output_file), 'velocity')
mp_sp_to_vtk (press, space_p, geometry, gnump, vtk_pts, sprintf ('%s_press', output_file), 'pressure')

