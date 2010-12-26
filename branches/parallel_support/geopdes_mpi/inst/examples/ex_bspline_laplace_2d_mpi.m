% EX_BSPLINE_LAPLACE_2D_MPI: Solve a 2d Laplace problem with a B-spline discretization. 
%
% Example to solve the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega = F((0,1)^2)
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% In order to make it work, a script file with the data must also be executed.
% The available data files for this problem are:
%
% - test_square
% - test_ring
% - test_ring_mixed_bc
% - test_plate_mixed_bc
%
% Copyright (C) 2009, 2010 Carlo de Falco
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Construct geometry structure
geometry  = geo_load (geo_name);

[knots, zeta] = kntrefine (geometry.nurbs.knots, n_sub, degree, regularity);
ser_zeta=zeta;

% initialize the parallel environment
Comm   = mpi_work_init();
Numm   = MPI_Comm_size(Comm);
Rank   = MPI_Comm_rank(Comm);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[zeta boundaries] = msh_partition (Numm, 'knots',knots,'degree',degree, 'index', Rank+1);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_2d_tensor_product (zeta, qn, qw, 'boundary', boundaries);
msh      = msh_push_forward_2d (msh, geometry);

[ser_qn,ser_qw] = msh_set_quad_nodes (ser_zeta, rule);
ser_msh  = msh_2d_tensor_product (ser_zeta, ser_qn, ser_qw);
ser_msh  = msh_push_forward_2d (ser_msh, geometry);

% Construct space structure
sp       = sp_bspline_2d_phys (knots, degree, msh);
ser_sp   = sp_bspline_2d_phys (knots, degree, ser_msh);

% Precompute the coefficients
x = squeeze (msh.geo_map(1,:,:));
y = squeeze (msh.geo_map(2,:,:));

ser_x = squeeze (ser_msh.geo_map(1,:,:));
ser_y = squeeze (ser_msh.geo_map(2,:,:));

epsilon = reshape (c_diff (x, y), msh.nqn, msh.nel);
fval    = reshape (f (x, y), msh.nqn, msh.nel) ;

ser_epsilon = reshape (c_diff (ser_x, ser_y), ser_msh.nqn, ser_msh.nel);
ser_fval    = reshape (f (ser_x, ser_y), ser_msh.nqn, ser_msh.nel);

% Assemble the matrices
stiff_mat = op_gradu_gradv (sp, sp, msh, epsilon);
rhs       = op_f_v (sp, msh, fval);

ser_stiff_mat = op_gradu_gradv (ser_sp, ser_sp, ser_msh, ser_epsilon);
ser_rhs       = op_f_v (ser_sp, ser_msh, ser_fval);

%norm(ser_rhs);
%norm(MPI_Allreduce(rhs,'mpi_sum',Comm)-ser_rhs)

% Apply Neumann boundary conditions
for iside = intersect( nmnn_sides, boundaries);
    x = squeeze (msh.boundary(iside).geo_map(1,:,:));
    y = squeeze (msh.boundary(iside).geo_map(2,:,:));
    gval = reshape (g (x, y, iside), msh.boundary(iside).nqn, msh.boundary(iside).nel);
  
    rhs(sp.boundary(iside).dofs) = rhs(sp.boundary(iside).dofs) + ...
        op_f_v (sp.boundary(iside), msh.boundary(iside), gval);
end

% Apply Dirichlet boundary conditions
u = zeros (sp.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj_mpi (sp, msh, h, drchlt_sides,Comm);
u(drchlt_dofs) = u_drchlt; 

int_dofs = setdiff (1:sp.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt;
rhs        = MPI_Allreduce(rhs, 'mpi_sum',Comm);

% parameters of the restarted parallel generalized minimmal residual solver
tol = 10^(-10);
guess = rhs(int_dofs);
max_iter = length(int_dofs);
restart  = min(50,length(int_dofs)); 

left_preconditioner  =  mpi_diagonal_preconditioner_fast( stiff_mat(int_dofs,int_dofs), Comm );

% Solve the linear system
u(int_dofs) = mpi_Rgmres(stiff_mat(int_dofs,int_dofs), rhs(int_dofs), guess, tol, max_iter, restart, Comm, {left_preconditioner,[]});
clear left_preconditioner;

% Postprocessing
if (exist ('uex', 'var'))
  error_l2 = sp_l2_error (sp, msh, u, uex);
  error_l2 = sqrt(MPI_Allreduce(error_l2^2,'mpi_sum',Comm));
  mpi_message(sprintf('L2 error: %g', error_l2),Comm);
  if (exist ('graduex', 'var'))
    error_h1 = sp_h1_error (sp, msh, u, uex, graduex);
    error_h1 = sqrt(MPI_Allreduce(error_h1^2,'mpi_sum',Comm));
    mpi_message(sprintf('L2 error: %g', error_h1),Comm);
  end  
end

% Serial post processing
if Rank ~= 0
  MPI_Finalize();
else
  fprintf ('The result is saved in the file %s \n \n', output_file);
  sp_to_vtk_2d (u, sp, geometry, vtk_pts, output_file, 'u');
  MPI_Finalize();
end





%!demo
%! test_square
%! ex_bspline_laplace_2d

%!demo
%! test_ring
%! ex_bspline_laplace_2d

%!demo
%! test_ring_mixed_bc
%! ex_bspline_laplace_2d

%!demo
%! test_plate_mixed_bc
%! ex_bspline_laplace_2d
