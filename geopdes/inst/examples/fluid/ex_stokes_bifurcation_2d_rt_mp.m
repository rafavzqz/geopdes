% EX_STOKES_BIFURCATION_RT_MP: data file for Stokes problem in a pipe with a bifurcation.

% 1) PHYSICAL DATA OF THE PROBLEM
problem_data  = struct ();

% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_bifurcation_mp.txt';

% Type of boundary conditions for each side of the domain
problem_data.drchlt_sides = 1:3;
problem_data.nmnn_sides = [];

% Physical parameters
problem_data.viscosity = @(x, y) ones (size (x));

% Force term
problem_data.f  = @(x, y) zeros ([2, size(x)]);

% Boundary terms
problem_data.h  = @test_stokes_bifurcation_mp_h_drchlt;

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
method_data = struct ();

method_data.element_name = 'rt';     % Element type for discretization

method_data.degree       = [3 3];  % Degree of the splines
method_data.regularity   = [2 2];  % Regularity of the splines
method_data.nsub         = [5 5];  % Number of subdivisions
method_data.nquad        = [5 5];  % Points for the Gaussian quadrature rule
method_data.Cpen = 10 * (method_data.degree(1)+1);

% 3) CALL TO THE SOLVER
[geometry, msh, space_v, vel, space_p, press] = mp_solve_stokes_div_conforming (problem_data, method_data);


% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW
output_file  = 'bifurcation_2d_rt_mp_deg3_reg2_sub5';
vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};

fprintf ('results being saved in: %s_vel.pvd and %s_press.pvd\n', output_file, output_file)
sp_to_vtk (vel, space_v, geometry, vtk_pts, sprintf ('%s_vel', output_file), {'velocity', 'divergence'}, {'value', 'divergence'})
sp_to_vtk (press, space_p, geometry, vtk_pts, sprintf ('%s_press', output_file), 'pressure')

%!test
%! problem_data  = struct ();
%! problem_data.geo_name = 'geo_bifurcation_mp.txt';
%! problem_data.drchlt_sides = 1:3;
%! problem_data.nmnn_sides = [];
%! problem_data.viscosity = @(x, y) ones (size (x));
%! problem_data.f  = @(x, y) zeros ([2, size(x)]);
%! problem_data.h  = @test_stokes_bifurcation_mp_h_drchlt;
%! method_data = struct ();
%! method_data.element_name = 'rt';     % Element type for discretization
%! method_data.degree       = [3 3];  % Degree of the splines
%! method_data.regularity   = [2 2];  % Regularity of the splines
%! method_data.nsub         = [5 5];  % Number of subdivisions
%! method_data.nquad        = [5 5];  % Points for the Gaussian quadrature rule
%! method_data.Cpen = 10 * (method_data.degree(1)+1);
%! [geometry, msh, space_v, vel, space_p, press] = mp_solve_stokes_div_conforming (problem_data, method_data);
%! assert (msh.nel, 100)
%! assert (space_v.ndof, 552)
%! assert (space_p.ndof, 256)
%! for iptc = 1:4
%!   div = sp_eval (vel(space_v.gnum{iptc}) .* space_v.dofs_ornt{iptc}', space_v.sp_patch{iptc}, geometry(iptc), [20 20], 'divergence');
%!   assert (max (abs (div(:))) < 1e-12)
%! end
