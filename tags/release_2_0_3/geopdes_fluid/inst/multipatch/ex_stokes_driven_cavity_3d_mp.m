% EX_STOKES_SQUARE_MP: data file for Stokes problem in the square defined with a multipatch geometry.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_2cubesa.txt';
% You can see how the patches should match in the files
%  geo_2cubesa{b,c,d,e,f,g,h}.txt
% In all the eight cases the result must be the same

% Type of boundary conditions for each side of the domain
problem_data.drchlt_sides = 1:6;
problem_data.nmnn_sides = [];

% Physical parameters
problem_data.viscosity = @(x, y, z) ones (size (x));

% Force term
problem_data.f  = @(x, y, z) zeros ([3, size(x)]);

% Boundary terms
problem_data.h  = @test_stokes_3d_symdrivcav_h_drchlt;

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.element_name = 'th';     % Element type for discretization
method_data.degree       = [2 2 2];  % Degree of the splines
method_data.regularity   = [1 1 1];  % Regularity of the splines
method_data.nsub         = [2 2 2];  % Number of subdivisions
method_data.nquad        = [4 4 4];  % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space_v, vel, gnum, space_p, press, gnump] = ...
                       mp_solve_stokes_3d (problem_data, method_data);
%[geometry, msh, space_v, vel, space_p, press] = ...
%                       solve_stokes_2d_bsplines (problem_data, method_data);

%  error_l2_p = sp_l2_error (space_p, msh, press, problem_data.pressex)
%  [error_h1, error_l2] = sp_h1_error (space_v, msh, ...
%     vel, problem_data.velex, problem_data.gradvelex)
%keyboard

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW
output_file  = 'DrivenCavity_3D_MP_Deg2_Reg1_Sub2';
vtk_pts = {linspace(0, 1, 10), linspace(0, 1, 10), linspace(0, 1, 10)};

fprintf ('results being saved in: %s_vel.pvd and %s_press.pvd\n', output_file, output_file)
mp_sp_to_vtk (vel, space_v, geometry, gnum, vtk_pts, sprintf ('%s_vel', output_file), 'velocity')
mp_sp_to_vtk (press, space_p, geometry, gnump, vtk_pts, sprintf ('%s_press', output_file), 'pressure')

