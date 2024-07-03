% PHYSICAL DATA OF THE PROBLEM
clear problem_data
problem_data.geo_name = 'geo_open_quasisphere_5p_ASG1.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];
problem_data.weak_drchlt_sides = [];

% Physical parameters
problem_data.c_diff  = @(x, y, z) ones(size(x));

% Source and boundary terms
problem_data.f = @(x, y, z) ones(size(x));
problem_data.g = @(x, y, z, ind) zeros(size(x));
problem_data.h = @(x, y, z, ind) zeros(size(x));

% DISCRETIZATION PARAMETERS
clear method_data
deg = 4;
method_data.degree      = [deg deg];     % Degree of the splines
method_data.regularity  = [deg-2 deg-2]; % Regularity of the splines
method_data.nsub        = [8 8];         % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nquad       = [deg+1 deg+1]; % Points for the Gaussian quadrature rule

[geometry, msh, space, u] = mp_solve_bilaplace_C1 (problem_data, method_data);
npts = [41 41];
sp_plot_solution (u, space, geometry); shading interp; axis equal
