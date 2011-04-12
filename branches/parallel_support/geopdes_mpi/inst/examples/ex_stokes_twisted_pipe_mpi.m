% EX_STOKES_TWISTED_PIPE: solve the Stokes problem on a twisted pipe using a bspline discretization.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_twisted_pipe.mat';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [1 2];
problem_data.drchlt_sides = [3 4 5 6];

% Physical parameters
problem_data.viscosity = @(x, y, z) ones (size (x));

% Force term
fx = @(x, y, z) ones(size(x));
fy = @(x, y, z) zeros(size(x));
fz = @(x, y, z) zeros(size(x));
problem_data.f  = @(x, y, z) cat(1, reshape (fx (x,y,z), [1, size(x)]), reshape (fy (x,y,z), [1, size(x)]), reshape (fz (x,y,z), [1, size(x)]));

% Boundary terms
problem_data.h  = @(x, y, z, iside) zeros ([3, size(x)]); %Dirichlet boundary condition
problem_data.g  = @(x, y, z, iside) zeros ([3, size(x)]); %Neumann boundary condition

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.element_name = 'TH';
method_data.degree       = [2 2 2]; % Degree of the splines
method_data.regularity   = [1 1 1]; % Regularity of the splines
method_data.nsub         = [2 2 2]; % Number of subdivisions
method_data.nquad        = [4 4 4]; % Points for the Gaussian quadrature rule

% Main solver
% left and right preconditioner
method_data.mpi_prec = @(matrix, comm) {mpi_diagonal_preconditioner_fast( matrix, comm ), [] }; 
% solver mpi_Rgmres (matrix, rhs, guess, tolerance, max_iterations, restart_interval, mpi_comm)
method_data.mpi_solv = @(matrix, rhs, prec, comm) mpi_Rgmres( matrix, rhs, rhs, 10^-10, size(matrix,1), min(50,size(matrix,1)) , comm, prec);

% Dirichlet condition solver
method_data.drchlt_prec = @(matrix, comm) mpi_diagonal_preconditioner_fast( matrix, comm );
method_data.drchlt_solv = @(matrix, rhs, prec, comm) mpi_Rgmres( matrix, rhs, rhs, 10^-10, size(matrix,1), min(50,size(matrix,1)) , comm, prec);

% 3)MPI Setup
mpi_comm  = mpi_work_init();  % Initialize MPI and choose a communication
			      % group for solving see MPI, and openmpi_ext
                              % documentation

% 4) CALL TO THE SOLVER
[geometry, msh, space_v, vel, space_p, press] = mpi_solve_stokes_3d (mpi_comm, problem_data, method_data);

% 5) POST-PROCESSING
% 5.1) EXPORT TO PARAVIEW


if MPI_Comm_rank(mpi_comm) == 0
  output_file = 'TwistedPipe_BSP_Deg2_Reg1_Sub2';
  fprintf ('The result is saved in the files %s \n and %s \n \n', ...
           [output_file '_vel'], [output_file '_press']);
  vtk_pts = {linspace(0, 1, 10), linspace(0, 1, 10), linspace(0, 1, 10)};
  sp_to_vtk_3d (press, space_p, geometry, vtk_pts, [output_file '_press'], 'press')
  sp_to_vtk_3d (vel,   space_v, geometry, vtk_pts, [output_file '_vel'  ], 'vel')
end

