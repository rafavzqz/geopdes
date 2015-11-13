% EX_STOKES_TWISTED_PIPE_RT_MP: data file for Stokes problem in the twisted pipe.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data

% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_twisted_pipe_mp.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides = [1 2];
problem_data.drchlt_sides = 3;

% Physical parameters
problem_data.viscosity = @(x, y, z) ones (size (x));

% Force term
fx = @(x, y, z) ones(size(x));
fy = @(x, y, z) zeros(size(x));
fz = @(x, y, z) zeros(size(x));
problem_data.f  = @(x, y, z) cat(1, reshape (fx (x,y,z), [1, size(x)]), reshape (fy (x,y,z), [1, size(x)]), reshape (fz (x,y,z), [1, size(x)]));

% Boundary terms
problem_data.h  = @(x, y, z, iside) zeros ([3, size(x)]); %Dirichlet boundary condition
problem_data.g  = @(x, y, z, iside) zeros ([3, size(x)]); %Neumann boundary condition

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.element_name = 'rt';     % Element type for discretization
method_data.degree       = [2 2 2];  % Degree of the splines
method_data.regularity   = [1 1 1];  % Regularity of the splines
method_data.nsub         = [2 2 2];  % Number of subdivisions
method_data.nquad        = [4 4 4];  % Points for the Gaussian quadrature rule
method_data.Cpen = 10 * (method_data.degree(1)+1);

% 3) CALL TO THE SOLVER
[geometry, msh, space_v, vel, space_p, press] = mp_solve_stokes_div_conforming (problem_data, method_data);


% 4) POST-PROCESSING
% EXPORT TO PARAVIEW
output_file  = 'TwistedPipe_RT';
vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20), linspace(0, 1, 20)};

fprintf ('results being saved in: %s_vel.pvd and %s_press.pvd\n', output_file, output_file)
sp_to_vtk (vel, space_v, geometry, vtk_pts, sprintf ('%s_vel', output_file), {'velocity', 'divergence'}, {'value', 'divergence'})
sp_to_vtk (press, space_p, geometry, vtk_pts, sprintf ('%s_press', output_file), 'pressure')


%!test
%! problem_data.geo_name = 'geo_twisted_pipe_mp.txt';
%! problem_data.nmnn_sides = [1 2];
%! problem_data.drchlt_sides = 3;
%! problem_data.viscosity = @(x, y, z) ones (size (x));
%! fx = @(x, y, z) ones(size(x));
%! fy = @(x, y, z) zeros(size(x));
%! fz = @(x, y, z) zeros(size(x));
%! problem_data.f  = @(x, y, z) cat(1, reshape (fx (x,y,z), [1, size(x)]), reshape (fy (x,y,z), [1, size(x)]), reshape (fz (x,y,z), [1, size(x)]));
%! problem_data.h  = @(x, y, z, iside) zeros ([3, size(x)]); %Dirichlet boundary condition
%! problem_data.g  = @(x, y, z, iside) zeros ([3, size(x)]); %Neumann boundary condition
%! method_data.element_name = 'rt';     % Element type for discretization
%! method_data.degree       = [2 2 2];  % Degree of the splines
%! method_data.regularity   = [1 1 1];  % Regularity of the splines
%! method_data.nsub         = [2 2 2];  % Number of subdivisions
%! method_data.nquad        = [4 4 4];  % Points for the Gaussian quadrature rule
%! method_data.Cpen = 10 * (method_data.degree(1)+1);
%! [geometry, msh, space_v, vel, space_p, press] = mp_solve_stokes_div_conforming (problem_data, method_data);
%! assert (msh.nel, 192)
%! assert (space_v.ndof, 3430)
%! assert (space_p.ndof, 1029)
%! for iptc = 1:3
%!   div = sp_eval (vel(space_v.gnum{iptc}) .* space_v.dofs_ornt{iptc}', space_v.sp_patch{iptc}, geometry(iptc), [20 20 20], 'divergence');
%!   assert ((max(abs(div(:))) - min(abs(div(:)))) < 1e-7)
%! end
