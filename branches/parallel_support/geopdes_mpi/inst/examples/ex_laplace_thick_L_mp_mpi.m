% EX_LAPLACE_THICK_L_MP: solve the Poisson problem in the multipatch thick L-shaped domain with a B-spline discretization.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
%  In both cases the result should be the same
problem_data.geo_name = 'geo_thickL_mp.txt';
%problem_data.geo_name = 'geo_thickL_mp_b.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [1 2];
problem_data.drchlt_sides = [3 4 5 6 7 8];

% Physical parameters
problem_data.c_diff  = @(x, y, z) ones(size(x));

% Source and boundary terms
problem_data.f = @(x, y, z) exp(x).*cos(z).*(-2*cos(x.*y).*y + sin(x.*y).*(y.^2+x.^2));
problem_data.g = @(x, y, z, ind) test_thick_Lshaped_mp_g_nmnn (x, y, z, ind);
problem_data.h = @(x, y, z, ind) exp (x) .* sin (x.*y) .* cos (z);

% Exact solution (optional)
problem_data.uex     = @(x, y, z) exp (x) .* sin (x.*y) .* cos (z);
problem_data.graduex = @(x, y, z) cat (1, ...
               reshape (exp (x) .* cos (z) .* (sin (x.*y) + y .* cos (x.*y)), [1, size(x)]), ...
               reshape (exp (x) .* x .* cos (x.*y) .* cos(z), [1, size(x)]), ...
               reshape (-exp (x) .* sin (x.*y) .* sin (z), [1, size(x)]));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [2 2 2];  % Degree of the splines
method_data.regularity = [1 1 1];  % Regularity of the splines
method_data.nsub       = [8 8 8];  % Number of subdivisions
method_data.nquad      = [3 3 3];  % Points for the Gaussian quadrature rule
% Main solver
% left and right preconditioner
method_data.mpi_prec = @(matrix, comm) {mpi_diagonal_preconditioner_fast( matrix, comm ), [] }; 
% solver mpi_Rgmres (matrix, rhs, guess, tolerance, max_iterations, restart_interval, mpi_comm)
method_data.mpi_solv = @(matrix, rhs, prec, comm) mpi_Rgmres( matrix, rhs, rhs, 10^-10, size(matrix,1), min(50,size(matrix,1)) , comm, prec);

% Dirichlet condition solver
method_data.drchlt_prec = @(matrix, comm) {mpi_diagonal_preconditioner_fast( matrix, comm ),[]};
method_data.drchlt_solv = @(matrix, rhs, prec, comm) mpi_Rgmres( matrix, rhs, rhs, 10^-10, size(matrix,1), min(50,size(matrix,1)) , comm, prec);

% 4) MPI Setup
mpi_comm  = mpi_work_init();  % Initialize MPI and choose a communication
			      % group for solving see MPI, and openmpi_ext
                              % documentation


% 5) CALL TO THE SOLVER
[geometry, msh, space, u, gnum] = ...
               mp_mpi_solve_laplace_3d (mpi_comm, problem_data, method_data);

% 6) POST-PROCESSING
% 6.1) EXPORT TO PARAVIEW
if MPI_Comm_rank(mpi_comm)==0
  output_file = 'thickL_mp_BSP_Deg2_Reg1_Sub8';
  vtk_pts = {linspace(0, 1, 10), linspace(0, 1, 10), linspace(0, 1, 10)};
  fprintf ('The result is saved in the file %s.pvd \n \n', output_file);
  mp_sp_to_vtk_3d (u, space, geometry, gnum, vtk_pts, output_file, 'u')
end
% 6.2) PLOT IN MATLAB. COMPARISON WITH THE EXACT SOLUTION
npatch = numel (geometry);
for iptc = 1:npatch
  [error_h1(iptc), error_l2(iptc)] = mpi_sp_h1_error (mpi_comm, space{iptc}, msh{iptc}, ...
     u(gnum{iptc}), problem_data.uex, problem_data.graduex);
end
error_l2 = sqrt (sum (error_l2 .* error_l2));
error_h1 = sqrt (sum (error_h1 .* error_h1));

mpi_message(sprintf('error L2 norm: %e\nerror H1 norm: %e\n', error_l2, error_h1), mpi_comm);

MPI_Finalize();

%!demo
%! ex_laplace_thick_L_mp_mpi