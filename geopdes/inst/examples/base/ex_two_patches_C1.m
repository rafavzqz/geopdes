% EX_LAPLACE_LSHAPED_MP: solve the Poisson problem in the multipatch L-shaped domain with a B-spline discretization.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
% nrb(1) = nrb4surf ([0 0], [-1 0], [0 1], [-1 1]);
% nrb(2) = nrb4surf ([0 0], [1 0], [0 1], [1 1]);

p1 = [-4 -1/2]; p2 = [0 0]; p3 = [-10/3 16/5]; p4 = [0 3];
nrb(1) = nrb4surf (p2, p1, p4, p3);

p1 = [0 0]; p2 = [8/3 -2/5]; p3 = [0 3]; p4 = [10/3 23/7];
nrb(2) = nrb4surf (p1, p2, p3, p4); problem_data.geo_id = 1;
problem_data.geo_name = nrb;

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [];
problem_data.weak_drchlt_sides = [1 2 3 4 5 6];

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

             
% problem_data.f = @(x, y) zeros (size(x));
% problem_data.h = @(x, y, ind) exp(x) .* sin (y);
% problem_data.uex     = @(x, y) exp (x) .* sin (y);
% problem_data.graduex = @(x, y) cat (1, ...
%                        reshape (exp(x).*sin(y), [1, size(x)]), ...
%                        reshape (exp(x).*cos(y), [1, size(x)]));
             
% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3];       % Degree of the splines
method_data.regularity = [1 1];       % Regularity of the splines
method_data.nsub       = [3 3];       % Number of subdivisions
method_data.nquad      = method_data.degree+1;       % Points for the Gaussian quadrature rule
method_data.Cpen = 10 * (min(method_data.degree) + 1);

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = mp_solve_laplace_C1 (problem_data, method_data);

% % 4) POST-PROCESSING
% % EXPORT TO PARAVIEW
% output_file = 'Lshaped_mp_BSP_Deg3_Reg2_Sub9';
% 
vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
% fprintf ('The result is saved in the file %s.pvd \n \n', output_file);
% sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')
% 
figure
sp_plot_solution (u, space, geometry, vtk_pts)

% COMPARISON WITH THE EXACT SOLUTION
[error_h1, error_l2] = sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex)

