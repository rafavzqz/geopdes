% EX_LAPLACE_CUBE_MP: solve the Poisson problem in the multipatch unit cube with a B-splines discretization.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_2cubesa.txt';
% You can see how the patches should match in the files
%  geo_2cubesa{b,c,d,e,f,g,h}.txt
% In all the eight cases the result must be the same

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4 5 6];

% Physical parameters
problem_data.c_diff  = @(x, y, z) ones(size(x));

% Source and boundary terms
problem_data.f = @(x, y, z) -exp (x + z) .* sin (y);
problem_data.h = @(x, y, z, ind) exp (x + z) .* sin (y); 

% Exact solution (optional)
problem_data.uex     = @(x, y, z) exp (x + z) .* sin (y);
problem_data.graduex = @(x, y, z) cat (1, ...
                          reshape (exp (x + z) .* sin (y), [1, size(x)]), ...
                          reshape (exp (x + z) .* cos (y), [1, size(x)]), ...
                          reshape (exp (x + z) .* sin (y), [1, size(x)]));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [2 2 2];  % Degree of the splines
method_data.regularity = [1 1 1];  % Regularity of the splines
method_data.nsub       = [3 3 3];  % Number of subdivisions
method_data.nquad      = [3 3 3];  % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = mp_solve_laplace (problem_data, method_data);

% 4) POST-PROCESSING
% EXPORT TO PARAVIEW
output_file = 'cube_mp_BSP_Deg2_Reg1_Sub3';

vtk_pts = {linspace(0, 1, 10), linspace(0, 1, 15), linspace(0, 1, 15)};
fprintf ('The result is saved in the file %s.pvd \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% COMPARISON WITH THE EXACT SOLUTION
[error_h1, error_l2] = sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex)

%!demo
%! ex_laplace_cube_mp

%!test
%! problem_data.geo_name = 'geo_2cubesa.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_diff  = @(x, y, z) ones(size(x));
%! problem_data.f = @(x, y, z) -exp (x + z) .* sin (y);
%! problem_data.h = @(x, y, z, ind) exp (x + z) .* sin (y); 
%! problem_data.uex     = @(x, y, z) exp (x + z) .* sin (y);
%! problem_data.graduex = @(x, y, z) cat (1, ...
%!                           reshape (exp (x + z) .* sin (y), [1, size(x)]), ...
%!                           reshape (exp (x + z) .* cos (y), [1, size(x)]), ...
%!                           reshape (exp (x + z) .* sin (y), [1, size(x)]));
%! method_data.degree     = [2 2 2];  % Degree of the splines
%! method_data.regularity = [1 1 1];  % Regularity of the splines
%! method_data.nsub       = [3 3 3];  % Number of subdivisions
%! method_data.nquad      = [3 3 3];  % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u, gnum] = mp_solve_laplace (problem_data, method_data);
%! for iptc = 1:numel(geometry)
%!   [error_h1(iptc), error_l2(iptc)] = sp_h1_error (space{iptc}, msh{iptc}, ...
%!      u(gnum{iptc}), problem_data.uex, problem_data.graduex);
%! end
%! error_l2 = sqrt (sum (error_l2 .* error_l2));
%! error_h1 = sqrt (sum (error_h1 .* error_h1));
%! assert (error_l2, 3.74776528960627e-04, 1e-16)
%! assert (error_h1, 0.0112329701133678, 1e-15)
%! assert (max ([gnum{:}]), numel (u))
%! assert (max ([gnum{:}]), 360)
%! assert (gnum{1}(space{1}.boundary(2).dofs), gnum{2}(space{2}.boundary(1).dofs))

%!test
%! problem_data.geo_name = 'geo_2cubesb.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_diff  = @(x, y, z) ones(size(x));
%! problem_data.f = @(x, y, z) -exp (x + z) .* sin (y);
%! problem_data.h = @(x, y, z, ind) exp (x + z) .* sin (y); 
%! problem_data.uex     = @(x, y, z) exp (x + z) .* sin (y);
%! problem_data.graduex = @(x, y, z) cat (1, ...
%!                           reshape (exp (x + z) .* sin (y), [1, size(x)]), ...
%!                           reshape (exp (x + z) .* cos (y), [1, size(x)]), ...
%!                           reshape (exp (x + z) .* sin (y), [1, size(x)]));
%! method_data.degree     = [2 2 2];  % Degree of the splines
%! method_data.regularity = [1 1 1];  % Regularity of the splines
%! method_data.nsub       = [3 3 3];  % Number of subdivisions
%! method_data.nquad      = [3 3 3];  % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u, gnum] = mp_solve_laplace (problem_data, method_data);
%! for iptc = 1:numel(geometry)
%!   [error_h1(iptc), error_l2(iptc)] = sp_h1_error (space{iptc}, msh{iptc}, ...
%!      u(gnum{iptc}), problem_data.uex, problem_data.graduex);
%! end
%! error_l2 = sqrt (sum (error_l2 .* error_l2));
%! error_h1 = sqrt (sum (error_h1 .* error_h1));
%! assert (error_l2, 3.74776528960627e-04, 1e-16)
%! assert (error_h1, 0.0112329701133678, 1e-15)
%! assert (max ([gnum{:}]), numel (u))
%! assert (max ([gnum{:}]), 360)
%! assert (gnum{1}(space{1}.boundary(2).dofs), gnum{2}(reshape (space{2}.boundary(1).dofs, space{2}.boundary(1).ndof_dir)'(:)))

%!test
%! problem_data.geo_name = 'geo_2cubesc.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_diff  = @(x, y, z) ones(size(x));
%! problem_data.f = @(x, y, z) -exp (x + z) .* sin (y);
%! problem_data.h = @(x, y, z, ind) exp (x + z) .* sin (y); 
%! problem_data.uex     = @(x, y, z) exp (x + z) .* sin (y);
%! problem_data.graduex = @(x, y, z) cat (1, ...
%!                           reshape (exp (x + z) .* sin (y), [1, size(x)]), ...
%!                           reshape (exp (x + z) .* cos (y), [1, size(x)]), ...
%!                           reshape (exp (x + z) .* sin (y), [1, size(x)]));
%! method_data.degree     = [2 2 2];  % Degree of the splines
%! method_data.regularity = [1 1 1];  % Regularity of the splines
%! method_data.nsub       = [3 3 3];  % Number of subdivisions
%! method_data.nquad      = [3 3 3];  % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u, gnum] = mp_solve_laplace (problem_data, method_data);
%! for iptc = 1:numel(geometry)
%!   [error_h1(iptc), error_l2(iptc)] = sp_h1_error (space{iptc}, msh{iptc}, ...
%!      u(gnum{iptc}), problem_data.uex, problem_data.graduex);
%! end
%! error_l2 = sqrt (sum (error_l2 .* error_l2));
%! error_h1 = sqrt (sum (error_h1 .* error_h1));
%! assert (error_l2, 3.74776528960627e-04, 1e-16)
%! assert (error_h1, 0.0112329701133678, 1e-15)
%! assert (max ([gnum{:}]), numel (u))
%! assert (max ([gnum{:}]), 360)
%! assert (gnum{1}(space{1}.boundary(2).dofs), gnum{2}(fliplr(reshape (space{2}.boundary(1).dofs, space{2}.boundary(1).ndof_dir))'(:)))

%!test
%! problem_data.geo_name = 'geo_2cubesd.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_diff  = @(x, y, z) ones(size(x));
%! problem_data.f = @(x, y, z) -exp (x + z) .* sin (y);
%! problem_data.h = @(x, y, z, ind) exp (x + z) .* sin (y); 
%! problem_data.uex     = @(x, y, z) exp (x + z) .* sin (y);
%! problem_data.graduex = @(x, y, z) cat (1, ...
%!                           reshape (exp (x + z) .* sin (y), [1, size(x)]), ...
%!                           reshape (exp (x + z) .* cos (y), [1, size(x)]), ...
%!                           reshape (exp (x + z) .* sin (y), [1, size(x)]));
%! method_data.degree     = [2 2 2];  % Degree of the splines
%! method_data.regularity = [1 1 1];  % Regularity of the splines
%! method_data.nsub       = [3 3 3];  % Number of subdivisions
%! method_data.nquad      = [3 3 3];  % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u, gnum] = mp_solve_laplace (problem_data, method_data);
%! for iptc = 1:numel(geometry)
%!   [error_h1(iptc), error_l2(iptc)] = sp_h1_error (space{iptc}, msh{iptc}, ...
%!      u(gnum{iptc}), problem_data.uex, problem_data.graduex);
%! end
%! error_l2 = sqrt (sum (error_l2 .* error_l2));
%! error_h1 = sqrt (sum (error_h1 .* error_h1));
%! assert (error_l2, 3.74776528960627e-04, 1e-16)
%! assert (error_h1, 0.0112329701133678, 1e-15)
%! assert (max ([gnum{:}]), numel (u))
%! assert (max ([gnum{:}]), 360)
%! assert (gnum{1}(space{1}.boundary(2).dofs), gnum{2}(flipud(reshape (space{2}.boundary(1).dofs, space{2}.boundary(1).ndof_dir))(:)))

%!test
%! problem_data.geo_name = 'geo_2cubese.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_diff  = @(x, y, z) ones(size(x));
%! problem_data.f = @(x, y, z) -exp (x + z) .* sin (y);
%! problem_data.h = @(x, y, z, ind) exp (x + z) .* sin (y); 
%! problem_data.uex     = @(x, y, z) exp (x + z) .* sin (y);
%! problem_data.graduex = @(x, y, z) cat (1, ...
%!                           reshape (exp (x + z) .* sin (y), [1, size(x)]), ...
%!                           reshape (exp (x + z) .* cos (y), [1, size(x)]), ...
%!                           reshape (exp (x + z) .* sin (y), [1, size(x)]));
%! method_data.degree     = [2 2 2];  % Degree of the splines
%! method_data.regularity = [1 1 1];  % Regularity of the splines
%! method_data.nsub       = [3 3 3];  % Number of subdivisions
%! method_data.nquad      = [3 3 3];  % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u, gnum] = mp_solve_laplace (problem_data, method_data);
%! for iptc = 1:numel(geometry)
%!   [error_h1(iptc), error_l2(iptc)] = sp_h1_error (space{iptc}, msh{iptc}, ...
%!      u(gnum{iptc}), problem_data.uex, problem_data.graduex);
%! end
%! error_l2 = sqrt (sum (error_l2 .* error_l2));
%! error_h1 = sqrt (sum (error_h1 .* error_h1));
%! assert (error_l2, 3.74776528960627e-04, 1e-16)
%! assert (error_h1, 0.0112329701133678, 1e-15)
%! assert (max ([gnum{:}]), numel (u))
%! assert (max ([gnum{:}]), 360)
%! assert (gnum{1}(space{1}.boundary(2).dofs), gnum{2}(fliplr(flipud(reshape (space{2}.boundary(1).dofs, space{2}.boundary(1).ndof_dir)))(:)))

%!test
%! problem_data.geo_name = 'geo_2cubesf.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_diff  = @(x, y, z) ones(size(x));
%! problem_data.f = @(x, y, z) -exp (x + z) .* sin (y);
%! problem_data.h = @(x, y, z, ind) exp (x + z) .* sin (y); 
%! problem_data.uex     = @(x, y, z) exp (x + z) .* sin (y);
%! problem_data.graduex = @(x, y, z) cat (1, ...
%!                           reshape (exp (x + z) .* sin (y), [1, size(x)]), ...
%!                           reshape (exp (x + z) .* cos (y), [1, size(x)]), ...
%!                           reshape (exp (x + z) .* sin (y), [1, size(x)]));
%! method_data.degree     = [2 2 2];  % Degree of the splines
%! method_data.regularity = [1 1 1];  % Regularity of the splines
%! method_data.nsub       = [3 3 3];  % Number of subdivisions
%! method_data.nquad      = [3 3 3];  % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u, gnum] = mp_solve_laplace (problem_data, method_data);
%! for iptc = 1:numel(geometry)
%!   [error_h1(iptc), error_l2(iptc)] = sp_h1_error (space{iptc}, msh{iptc}, ...
%!      u(gnum{iptc}), problem_data.uex, problem_data.graduex);
%! end
%! error_l2 = sqrt (sum (error_l2 .* error_l2));
%! error_h1 = sqrt (sum (error_h1 .* error_h1));
%! assert (error_l2, 3.74776528960627e-04, 1e-16)
%! assert (error_h1, 0.0112329701133678, 1e-15)
%! assert (max ([gnum{:}]), numel (u))
%! assert (max ([gnum{:}]), 360)
%! assert (gnum{1}(space{1}.boundary(2).dofs), gnum{2}(fliplr(flipud(reshape (space{2}.boundary(1).dofs, space{2}.boundary(1).ndof_dir)))'(:)))

%!test
%! problem_data.geo_name = 'geo_2cubesg.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_diff  = @(x, y, z) ones(size(x));
%! problem_data.f = @(x, y, z) -exp (x + z) .* sin (y);
%! problem_data.h = @(x, y, z, ind) exp (x + z) .* sin (y); 
%! problem_data.uex     = @(x, y, z) exp (x + z) .* sin (y);
%! problem_data.graduex = @(x, y, z) cat (1, ...
%!                           reshape (exp (x + z) .* sin (y), [1, size(x)]), ...
%!                           reshape (exp (x + z) .* cos (y), [1, size(x)]), ...
%!                           reshape (exp (x + z) .* sin (y), [1, size(x)]));
%! method_data.degree     = [2 2 2];  % Degree of the splines
%! method_data.regularity = [1 1 1];  % Regularity of the splines
%! method_data.nsub       = [3 3 3];  % Number of subdivisions
%! method_data.nquad      = [3 3 3];  % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u, gnum] = mp_solve_laplace (problem_data, method_data);
%! for iptc = 1:numel(geometry)
%!   [error_h1(iptc), error_l2(iptc)] = sp_h1_error (space{iptc}, msh{iptc}, ...
%!      u(gnum{iptc}), problem_data.uex, problem_data.graduex);
%! end
%! error_l2 = sqrt (sum (error_l2 .* error_l2));
%! error_h1 = sqrt (sum (error_h1 .* error_h1));
%! assert (error_l2, 3.74776528960627e-04, 1e-16)
%! assert (error_h1, 0.0112329701133678, 1e-15)
%! assert (max ([gnum{:}]), numel (u))
%! assert (max ([gnum{:}]), 360)
%! assert (gnum{1}(space{1}.boundary(2).dofs), gnum{2}(flipud(reshape (space{2}.boundary(1).dofs, space{2}.boundary(1).ndof_dir))'(:)))

%!test
%! problem_data.geo_name = 'geo_2cubesh.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_diff  = @(x, y, z) ones(size(x));
%! problem_data.f = @(x, y, z) -exp (x + z) .* sin (y);
%! problem_data.h = @(x, y, z, ind) exp (x + z) .* sin (y); 
%! problem_data.uex     = @(x, y, z) exp (x + z) .* sin (y);
%! problem_data.graduex = @(x, y, z) cat (1, ...
%!                           reshape (exp (x + z) .* sin (y), [1, size(x)]), ...
%!                           reshape (exp (x + z) .* cos (y), [1, size(x)]), ...
%!                           reshape (exp (x + z) .* sin (y), [1, size(x)]));
%! method_data.degree     = [2 2 2];  % Degree of the splines
%! method_data.regularity = [1 1 1];  % Regularity of the splines
%! method_data.nsub       = [3 3 3];  % Number of subdivisions
%! method_data.nquad      = [3 3 3];  % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u, gnum] = mp_solve_laplace (problem_data, method_data);
%! for iptc = 1:numel(geometry)
%!   [error_h1(iptc), error_l2(iptc)] = sp_h1_error (space{iptc}, msh{iptc}, ...
%!      u(gnum{iptc}), problem_data.uex, problem_data.graduex);
%! end
%! error_l2 = sqrt (sum (error_l2 .* error_l2));
%! error_h1 = sqrt (sum (error_h1 .* error_h1));
%! assert (error_l2, 3.74776528960627e-04, 1e-16)
%! assert (error_h1, 0.0112329701133678, 1e-15)
%! assert (max ([gnum{:}]), numel (u))
%! assert (max ([gnum{:}]), 360)
%! assert (gnum{1}(space{1}.boundary(2).dofs), gnum{2}(fliplr(reshape (space{2}.boundary(1).dofs, space{2}.boundary(1).ndof_dir))(:)))
