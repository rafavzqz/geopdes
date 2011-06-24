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
[geometry, msh, space, u, gnum] = ...
               mp_solve_laplace_3d (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW
output_file = 'cube_mp_BSP_Deg2_Reg1_Sub3';

vtk_pts = {linspace(0, 1, 10), linspace(0, 1, 15), linspace(0, 1, 15)};
fprintf ('The result is saved in the file %s.pvd \n \n', output_file);
mp_sp_to_vtk (u, space, geometry, gnum, vtk_pts, output_file, 'u')

% 4.2) PLOT IN MATLAB. COMPARISON WITH THE EXACT SOLUTION
npatch = numel (geometry);
for iptc = 1:npatch
  [error_h1(iptc), error_l2(iptc)] = sp_h1_error (space{iptc}, msh{iptc}, ...
     u(gnum{iptc}), problem_data.uex, problem_data.graduex);
end
error_l2 = sqrt (sum (error_l2 .* error_l2))
error_h1 = sqrt (sum (error_h1 .* error_h1))

%!demo
%! ex_laplace_cube_mp
