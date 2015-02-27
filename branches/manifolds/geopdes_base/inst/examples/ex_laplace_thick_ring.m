% EX_LAPLACE_THICK_RING: solve the Poisson problem in a thick ring with a B-spline discretization (non-isoparametric approach).

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_thick_ring.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [4 5 6];
problem_data.drchlt_sides = [1 2 3];

% Physical parameters
problem_data.c_diff  = @(x, y, z) ones(size(x));

% Source and boundary terms
problem_data.f = ...
     @(x, y, z) exp(x).*cos(z).*(-2*cos(x.*y).*y + sin(x.*y).*(y.^2 + x.^2));
problem_data.g = @test_thick_ring_g_nmnn;
problem_data.h = @(x, y, z, ind) exp (x) .* sin (x.*y) .* cos (z);

% Exact solution (optional)
problem_data.uex     = @(x, y, z) exp (x) .* sin (x.*y) .* cos (z);
problem_data.graduex = @(x, y, z) cat (1, ...
             reshape (exp (x) .* cos (z) .* (sin (x.*y) + y .* cos (x.*y)), [1, size(x)]), ...
             reshape (exp (x) .* x .* cos (x.*y) .* cos(z), [1, size(x)]), ...
             reshape (-exp (x) .* sin (x.*y) .* sin (z), [1, size(x)]));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [2 2 2];       % Degree of the splines
method_data.regularity = [1 1 1];       % Regularity of the splines
method_data.nsub       = [4 4 4];       % Number of subdivisions
method_data.nquad      = [3 3 3];       % Points for Gaussian quadrature rule

% 3) CALL TO THE SOLVER

[geometry, msh, space, u] = solve_laplace_3d (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW

output_file = 'Thick_ring_BSP_Deg2_Reg1_Sub4';

vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20), linspace(0, 1, 20)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% 4.2) COMPARISON WITH THE EXACT SOLUTION

[error_h1, error_l2] = ...
           sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex)

%!demo
%! ex_laplace_thick_ring
