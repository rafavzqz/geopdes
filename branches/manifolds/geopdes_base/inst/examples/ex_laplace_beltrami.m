% EX_LAPLACE_RING: solve the Poisson problem in one quarter of a ring, discretized with B-splines (non-isoparametric approach).

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS structure from the NURBS toolbox
problem_data.geo_name = 'geo_roof.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];

% Physical parameters
problem_data.c_diff  = @(x, y, z) ones(size(x));

% Source and boundary terms
alp = 5; bet = 1; L = 1;
g1 = @(x,y) (1-cos(atan(x./y))).*(1-sin(atan(x./y)));
g2 = @(x,y) cos(atan(x./y)) + sin(atan(x./y)) - 4*cos(atan(x./y)).*sin(atan(x./y));
gz = @(z) sin(pi*alp/L*z);
dg1 = @(x,y) sin(atan(y./x)).*(1-sin(atan(y./x))) - cos(atan(y./x)).*(1-cos(atan(y./x)));

problem_data.f = @(x, y, z) bet*((alp*pi/L)^2 * g1(x,y) - g2(x,y)).*gz(z);
problem_data.h = @(x, y, z, ind) bet * g1(x,y).*gz(z);

% Exact solution (optional)
problem_data.uex     = @(x, y, z) bet * g1(x,y).*gz(z);
problem_data.graduex = @(x, y, z) cat (1, ...
                reshape (-1./sqrt(x.^2+y.^2).*sin(atan(y./x))*bet.*gz(z).*dg1(x,y), [1, size(x)]), ...
                reshape ( 1./sqrt(x.^2+y.^2).*cos(atan(y./x))*bet.*gz(z).*dg1(x,y), [1, size(x)]), ...
                reshape ( alp*pi/L * cos(alp*pi/L*z) .* bet .* g1(x,y), [1, size(x)]));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3];       % Degree of the splines
method_data.regularity = [2 2];       % Regularity of the splines
method_data.nsub       = [9 9];       % Number of subdivisions
method_data.nquad      = [4 4];       % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_laplace (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW

output_file = 'Laplace-Beltrami.vts';

vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% 4.2) PLOT IN MATLAB. COMPARISON WITH THE EXACT SOLUTION

[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y, Z]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)), squeeze(F(3,:,:)));
subplot (1,2,1)
surf (X, Y, Z, eu)
title ('Numerical solution'), axis tight
subplot (1,2,2)
surf (X, Y, Z, problem_data.uex (X,Y,Z))
title ('Exact solution'), axis tight

% Display errors of the computed solution in the L2 and H1 norm
[error_h1, error_l2] = ...
           sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex)

%!demo
%! ex_laplace_beltrami

%!test
%! problem_data.geo_name = 'geo_roof.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4];
%! problem_data.c_diff  = @(x, y, z) ones(size(x));
%! alp = 5; bet = 1; L = 1;
%! g1 = @(x,y) (1-cos(atan(x./y))).*(1-sin(atan(x./y)));
%! g2 = @(x,y) cos(atan(x./y)) + sin(atan(x./y)) - 4*cos(atan(x./y)).*sin(atan(x./y));
%! gz = @(z) sin(pi*alp/L*z);
%! dg1 = @(x,y) sin(atan(y./x)).*(1-sin(atan(y./x))) - cos(atan(y./x)).*(1-cos(atan(y./x)));
%! problem_data.f = @(x, y, z) bet*((alp*pi/L)^2 * g1(x,y) - g2(x,y)).*gz(z);
%! problem_data.h = @(x, y, z, ind) bet * g1(x,y).*gz(z);
%! problem_data.uex     = @(x, y, z) bet * g1(x,y).*gz(z);
%! problem_data.graduex = @(x, y, z) cat (1, ...
%!                 reshape (-1./sqrt(x.^2+y.^2).*sin(atan(y./x))*bet.*gz(z).*dg1(x,y), [1, size(x)]), ...
%!                 reshape ( 1./sqrt(x.^2+y.^2).*cos(atan(y./x))*bet.*gz(z).*dg1(x,y), [1, size(x)]), ...
%!                 reshape ( alp*pi/L * cos(alp*pi/L*z) .* bet .* g1(x,y), [1, size(x)]));
%! method_data.degree     = [3 3];       % Degree of the splines
%! method_data.regularity = [2 2];       % Regularity of the splines
%! method_data.nsub       = [9 9];       % Number of subdivisions
%! method_data.nquad      = [4 4];       % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u] = solve_laplace (problem_data, method_data);
%! [error_h1, error_l2] = ...
%!            sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex);
%! assert (msh.nel, 81)
%! assert (space.ndof, 144)
%! assert (error_h1, 0.0407258650845198, 1e-15)
%! assert (error_l2, 0.00101798638120372, 1e-16)