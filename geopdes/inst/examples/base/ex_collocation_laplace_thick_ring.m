% EX_COLLOCATION_LAPLACE_THICK_RING: solve the Laplace problem with a NURBS discretization by isogeometric collocation

% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_thick_ring.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4 5 6];

% Physical parameters
problem_data.c_diff  = @(x, y, z) ones(size(x)); % Only used in Galerkin, not in collocation

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

% Discretization parameters (p and h)
method_data.degree     = [4 4 4]; % Degree of the splines in each direction
method_data.regularity = [3 3 3]; % Regularity of the splines, should be at least C^1.
method_data.nsub       = [4 4 4]; % Divide each subinterval of the original knot span in nsub subintervals
method_data.pts_case   = 1;       % Collocation points. 1: Greville abscissae, 2: clustered superconvergent points

% Solve with collocation method
[geometry, msh_coll, space_coll, u_coll] = solve_laplace_collocation (problem_data, method_data);

% Solve with Galerkin, for comparison
method_data.nquad = method_data.degree + 1;
[~, msh_gal, space_gal, u_gal] = solve_laplace_iso (problem_data, method_data);

% Plot of solution
figure;
subplot (1,2,1)
sp_plot_solution (u_coll, space_coll, geometry)
title ('Collocation solution'), axis tight
subplot (1,2,2)
sp_plot_solution (u_gal, space_gal, geometry)
title ('Galerkin solution'), axis tight


% Compute errors of collocation and Galerkin. Since it involves quadrature,
%  it is computed using the mesh and space objects from the Galerkin method.
disp ('Error in H1 and L2 norms, for isogeometric collocation')
[error_h1_coll, error_l2_coll] = sp_h1_error (space_gal, msh_gal, u_coll, problem_data.uex, problem_data.graduex);
disp([error_l2_coll error_h1_coll])

disp ('Error in H1 and L2 norms, for isogeometric Galerkin')
[error_h1_gal, error_l2_gal] = sp_h1_error (space_gal, msh_gal, u_gal, problem_data.uex, problem_data.graduex);
disp ([error_l2_gal error_h1_gal])
