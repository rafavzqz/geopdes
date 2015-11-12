% EX_LAPLACE_LSHAPED_MP: solve the Poisson problem in the multipatch L-shaped domain with a B-spline discretization.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_Lshaped_mp.txt';
%problem_data.geo_name = 'geo_Lshaped_mp_b.txt';
%problem_data.geo_name = 'geo_Lshaped_mp_c.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [3 4 5 6];
problem_data.drchlt_sides = [1 2];

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));

% Source and boundary terms
problem_data.f = @(x, y) exp(x).*((x.^2 + y.^2 - 1).*sin(x.*y) - 2*y.*cos(x.*y));
problem_data.g = @(x, y, ind) test_Lshaped_mp_g_nmnn (x, y, ind);
problem_data.h = @(x, y, ind) exp(x) .* sin (x.*y);

% Exact solution (optional)
problem_data.uex     = @(x, y) exp(x) .* sin (x.*y);
problem_data.graduex = @(x, y) cat (1, ...
               reshape (exp(x).*(sin(x.*y) + y.*cos(x.*y)), [1, size(x)]), ...
               reshape (exp(x).*x.*cos(x.*y), [1, size(x)]));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3];       % Degree of the splines
method_data.regularity = [2 2];       % Regularity of the splines
method_data.nsub       = [9 9];       % Number of subdivisions
method_data.nquad      = [4 4];       % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = mp_solve_laplace (problem_data, method_data);

% 4) POST-PROCESSING
% EXPORT TO PARAVIEW
output_file = 'Lshaped_mp_BSP_Deg3_Reg2_Sub9';

vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
fprintf ('The result is saved in the file %s.pvd \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% COMPARISON WITH THE EXACT SOLUTION
[error_h1, error_l2] = sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex)

%!demo
%! ex_laplace_Lshaped_mp

%!test
%! problem_data.geo_name = 'geo_Lshaped_mp.txt';
%! problem_data.nmnn_sides   = [3 4 5 6];
%! problem_data.drchlt_sides = [1 2];
%! problem_data.c_diff  = @(x, y) ones(size(x));
%! problem_data.f = @(x, y) exp(x).*((x.^2 + y.^2 - 1).*sin(x.*y) - 2*y.*cos(x.*y));
%! problem_data.g = @(x, y, ind) test_Lshaped_mp_g_nmnn (x, y, ind);
%! problem_data.h = @(x, y, ind) exp(x) .* sin (x.*y);
%! problem_data.uex     = @(x, y) exp(x) .* sin (x.*y);
%! problem_data.graduex = @(x, y) cat (1, ...
%!                reshape ( exp(x).*(sin(x.*y) + y.*cos(x.*y)), [1, size(x)]), ...
%!                reshape (exp(x).*x.*cos(x.*y), [1, size(x)]));
%! method_data.degree     = [3 3];       % Degree of the splines
%! method_data.regularity = [2 2];       % Regularity of the splines
%! method_data.nsub       = [9 9];       % Number of subdivisions
%! method_data.nquad      = [4 4];       % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u] = mp_solve_laplace (problem_data, method_data);
%! [error_h1, error_l2] = sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex);
%! assert (error_l2, 3.00642608570282e-07, 1e-16)
%! assert (error_h1, 1.77524941085757e-05, 1e-16)
%! assert (space.ndof, 408)
%! for iptc = 1:space.boundary.npatch
%!   patch = msh.boundary.patch_numbers(iptc);
%!   side = msh.boundary.side_numbers(iptc);
%!   assert (space.gnum{patch}(space.sp_patch{patch}.boundary(side).dofs)(:), space.boundary.dofs(space.boundary.gnum{iptc})(:))
%! end
%!
%! problem_data.geo_name = 'geo_Lshaped_mp_b.txt';
%! [geometry, msh, space, u] = mp_solve_laplace (problem_data, method_data);
%! [error_h1, error_l2] = sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex);
%! assert (error_l2, 3.00642608570282e-07, 1e-16)
%! assert (error_h1, 1.77524941085757e-05, 1e-16)
%! assert (space.ndof, 408)
%! for iptc = 1:space.boundary.npatch
%!   patch = msh.boundary.patch_numbers(iptc);
%!   side = msh.boundary.side_numbers(iptc);
%!   assert (space.gnum{patch}(space.sp_patch{patch}.boundary(side).dofs)(:), space.boundary.dofs(space.boundary.gnum{iptc})(:))
%! end
%!
%! problem_data.geo_name = 'geo_Lshaped_mp_c.txt';
%! [geometry, msh, space, u] = mp_solve_laplace (problem_data, method_data);
%! [error_h1, error_l2] = sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex);
%! assert (error_l2, 3.00642608570282e-07, 1e-16)
%! assert (error_h1, 1.77524941085757e-05, 1e-16)
%! assert (space.ndof, 408)
%! for iptc = 1:space.boundary.npatch
%!   patch = msh.boundary.patch_numbers(iptc);
%!   side = msh.boundary.side_numbers(iptc);
%!   assert (space.gnum{patch}(space.sp_patch{patch}.boundary(side).dofs)(:), space.boundary.dofs(space.boundary.gnum{iptc})(:))
%! end
