% EX_STOKES_DRIVEN_CAVITY_TH: solve the Stokes problem in the driven cavity with generalized Taylor-Hood elements.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined from the aspect ratio using the NURBS toolbox
%aspect_ratio=[ 2 3];
%temp_geo = geo_load('geo_cube.txt');
%problem_data.geo_name = nrbtform (temp_geo.nurbs, vecscale ([1 aspect_ratio(1) aspect_ratio(2)]));
problem_data.geo_name = 'geo_cube.txt';

% Type of boundary conditions for each side of the domain
problem_data.drchlt_sides = 1:6;
problem_data.nmnn_sides = [];

% Physical parameters
problem_data.viscosity = @(x, y, z) ones (size (x));

% Force and boundary terms
problem_data.f  = @(x, y, z) zeros ([3, size(x)]);
problem_data.h  = @test_stokes_3d_symdrivcav_h_drchlt;

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.element_name = 'sg';        % Element type for discretization
method_data.degree       = [ 2  2  2];  % Degree of the splines
method_data.regularity   = [ 1  1  1];  % Regularity of the splines
method_data.nsub         = [ 2  2  2];  % Number of subdivisions
method_data.nquad        = [ 4  4  4];  % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space_v, vel, space_p, press] = ...
                       solve_stokes_3d_bsplines (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW
output_file = 'Driven_cavity_3d_SG_Deg3_Reg2_Sub2'

vtk_pts = {linspace(0, 1, 10)', linspace(0, 1, 10)', linspace(0, 1, 10)'};
sp_to_vtk_3d (press, space_p, geometry, vtk_pts, [output_file '_press'], 'press')
sp_to_vtk_3d (vel,   space_v, geometry, vtk_pts, [output_file '_vel'  ], 'vel')

