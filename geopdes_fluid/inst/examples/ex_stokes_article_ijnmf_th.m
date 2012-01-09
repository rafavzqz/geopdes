% EX_STOKES_ARTICLE_IJNMF_TH: solve the Stokes problem in the driven cavity with generalized Taylor-Hood elements.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined from the aspect ratio using the NURBS toolbox
aspect_ratio = 0.9;
nrb_square = nrb4surf ([0 0], [1 0], [0 1], [1 1]);
problem_data.geo_name = nrbtform (nrb_square, vecscale ([1 aspect_ratio]));

% Type of boundary conditions for each side of the domain
problem_data.drchlt_sides = 1:4;
problem_data.nmnn_sides = [];

% Physical parameters
problem_data.viscosity = @(x, y) ones (size (x));

% Force and boundary terms
problem_data.f  = @(x, y) zeros ([2, size(x)]);
problem_data.h  = @test_stokes_symdrivcav_h_drchlt;

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.element_name = 'th';   % Element type for discretization
method_data.degree       = [ 3  3];  % Degree of the splines
method_data.regularity   = [ 2  2];  % Regularity of the splines
method_data.nsub         = [10 10];  % Number of subdivisions
method_data.nquad        = [ 5  5];  % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space_v, vel, space_p, press] = ...
                       solve_stokes_2d (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW
output_file = 'Driven_cavity_TH_Deg3_Reg2_Sub10';

fprintf ('The result is saved in the files %s \n and %s \n \n', ...
           [output_file '_vel'], [output_file '_press']);
vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
sp_to_vtk (press, space_p, geometry, vtk_pts, [output_file '_press'], 'press')
sp_to_vtk (vel,   space_v, geometry, vtk_pts, [output_file '_vel'  ], 'vel')

% 4.2) PLOT IN MATLAB
[eu, F] = sp_eval (vel, space_v, geometry, vtk_pts);
[X,  Y] = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
figure()
quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
axis equal
title('Computed solution')

[div, F] = sp_eval (vel, space_v, geometry, vtk_pts, 'divergence');
figure()
surf (X, Y, div)
view(2)
axis equal
title('Computed divergence')

%!demo
%! ex_stokes_article_ijnmf_th
