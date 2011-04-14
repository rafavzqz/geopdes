% EX_LAPLACE_CUBE: solve the Poisson problem in the unit cube with a B-spline discretization.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_cube.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [4 5 6];
problem_data.drchlt_sides = [1 2 3];

% Physical parameters
problem_data.c_diff  = @(x, y, z) ones(size(x));

% Source and boundary terms
problem_data.f = @(x, y, z) -exp(x+z).*sin(y);
problem_data.g = @test_cube_g_nmnn;
problem_data.h = @(x, y, z, ind) exp (x + z) .* sin (y);

% Exact solution (optional)
problem_data.uex     = @(x, y, z) exp (x + z) .* sin (y);
problem_data.graduex = @(x, y, z) cat (1, ...
                          reshape (exp (x + z) .* sin (y), [1, size(x)]), ...
                          reshape (exp (x + z) .* cos (y), [1, size(x)]), ...
                          reshape (exp (x + z) .* sin (y), [1, size(x)]));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree      = [3 3 3];       % Degree of the splines
method_data.regularity  = [2 2 2];       % Regularity of the splines
method_data.nsub        = [5 5 5];       % Number of subdivisions
method_data.nquad       = [4 4 4];       % Points for Gaussian quadrature rule

% Main solver
% left and right preconditioner
method_data.mpi_prec = @(matrix, comm) {mpi_diagonal_preconditioner_fast( matrix, comm ), [] }; 
% solver mpi_Rgmres (matrix, rhs, guess, tolerance, max_iterations, restart_interval, mpi_comm)
method_data.mpi_solv = @(matrix, rhs, prec, comm) mpi_Rgmres( matrix, rhs, rhs, 10^-10, size(matrix,1), min(50,size(matrix,1)) , comm, prec);

% Dirichlet condition solver
method_data.drchlt_prec = @(matrix, comm) {mpi_diagonal_preconditioner_fast( matrix, comm ), [] }; 
method_data.drchlt_solv = @(matrix, rhs, prec, comm) mpi_Rgmres( matrix, rhs, rhs, 10^-10, size(matrix,1), min(50,size(matrix,1)) , comm, prec);

% 3) MPI Setup
mpi_comm  = mpi_work_init();  % Initialize MPI and choose a communication
			      % group for solving see MPI, and openmpi_ext
                              % documentation

% 4) CALL TO THE SOLVER
[geometry, msh, space, u] = mpi_solve_laplace_3d (mpi_comm, problem_data, method_data);

% 5) POST-PROCESSING
% 5.1) EXPORT TO PARAVIEW

if MPI_Comm_rank(mpi_comm) == 0
  output_file = 'Cube_BSP_Deg3_Reg2_Sub2';
  vtk_pts = {linspace(0, 1, 10), linspace(0, 1, 10), linspace(0, 1, 10)};
  fprintf ('The result is saved in the file %s \n \n', output_file);
  sp_to_vtk_3d (u, space, geometry, vtk_pts, output_file, 'u')
end

% 5.2) COMPARISON WITH THE EXACT SOLUTION

[error_h1, error_l2] = ...
           mpi_sp_h1_error (mpi_comm, space, msh, u, problem_data.uex, problem_data.graduex);
mpi_message(sprintf('error L2 norm: %e\nerror H1 norm: %e\n', error_l2, error_h1), mpi_comm);

MPI_Finalize();

%!demo
%! ex_laplace_cube_mpi
