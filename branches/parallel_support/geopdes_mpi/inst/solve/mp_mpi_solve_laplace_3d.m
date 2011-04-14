% MP_SOLVE_LAPLACE_3D: solve the Laplacian problem in a multipatch geometry.
%
% Example to solve the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% where the domain \Omega is formed by several patches of the form F((0,1)^3).
%
% USAGE:
%
%  [geometry, msh, space, u, gnum] = 
%          mp_mpi_solve_laplace_3d (mpi_comm, problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_diff:       diffusion coefficient (epsilon in the equation)
%    - f:            source term
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%    - mpi_prec    preconditioner for the iterative method
%    - mpi_solv    solver for the system
%    - drchlt_prec as above but used for solving boundary conditions
%    - drchlt_solv
%
% OUTPUT:
%
%  geometry: array of geometry structures (see geo_load)
%  msh:      array of mesh structures (see msh_push_forward_2d)
%  space:    array of space structures (see sp_bspline_2d_phys)
%  u:        the computed degrees of freedom
%  gnum:     global numbering of the degrees of freedom
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2010, 2011 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 2 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [geometry, msh, sp, u, gnum] = ...
              mp_mpi_solve_laplace_3d (mpi_comm, problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} sprintf('= problem_data.(data_names{iopt});')]);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} sprintf('= method_data.(data_names{iopt});')]);
end

% Construct geometry structure, and information for interfaces and boundaries
[geometry, boundaries, interfaces] = mp_geo_load (geo_name);
npatch = numel (geometry);

msh = cell (1, npatch); 
sp  = cell (1, npatch);

for iptc = 1:npatch
% Define the refined mesh, with tensor product structure
  [knots{iptc}, zeta{iptc}] = ...
         kntrefine (geometry(iptc).nurbs.knots, nsub-1, degree, regularity);
end

% Identify ourselves and partition work
mpi_size   = MPI_Comm_size(mpi_comm);
mpi_rank   = MPI_Comm_rank(mpi_comm);
[mzeta boundary_list ptc] = mp_msh_partition (mpi_size, knots, degree, mpi_rank+1);

% Compute the quadrature rule
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (mzeta, rule);
msh{ptc} = msh_3d_tensor_product (mzeta, qn, qw, 'boundary', boundary_list);
msh{ptc} = msh_push_forward_3d   (msh{ptc}, geometry(ptc));
% Evaluate the discrete space basis functions in the quadrature points
sp{ptc} = sp_bspline_3d_phys (knots{ptc}, degree, msh{ptc});

rule     = msh_gauss_nodes ([0 0 0]);
for iptc = setdiff(1:npatch, ptc)
  % Collect the information needed to assign the global dofs numbers
  [qn, qw] = msh_set_quad_nodes (zeta{iptc}, rule);
  msh{iptc} = msh_3d_tensor_product (zeta{iptc}, qn, qw);
  msh{iptc} = msh_push_forward_3d (msh{iptc}, geometry(iptc));
  sp{iptc} = sp_bspline_3d_phys (knots{iptc}, degree, msh{iptc});
end

% Create a correspondence between patches on the interfaces
[gnum, ndof] = mp_interface_3d (interfaces, sp);

% Compute and assemble the matrices 
stiff_mat = spalloc (ndof, ndof, ndof);
rhs = zeros (ndof, 1);

for iptc = ptc

  x = squeeze (msh{iptc}.geo_map (1,:,:));
  y = squeeze (msh{iptc}.geo_map (2,:,:));
  z = squeeze (msh{iptc}.geo_map (3,:,:));

  epsilon = reshape (c_diff (x, y, z), msh{iptc}.nqn, msh{iptc}.nel);
  fval    = reshape (f(x, y, z), msh{iptc}.nqn, msh{iptc}.nel);

  A_loc   = op_gradu_gradv (sp{iptc}, sp{iptc}, msh{iptc}, epsilon);
  rhs_loc = op_f_v (sp{iptc}, msh{iptc}, fval);

  stiff_mat(gnum{iptc},gnum{iptc}) = stiff_mat(gnum{iptc},gnum{iptc}) + A_loc;
  rhs(gnum{iptc}) = rhs(gnum{iptc}) + rhs_loc;
end

for iref = nmnn_sides
  for bnd_side = 1:boundaries(iref).nsides
    iptc  = boundaries(iref).patches(bnd_side);
    if iptc ~= ptc
      continue
    end
    iside = boundaries(iref).faces(bnd_side);
    if ~ismember(iside, msh{iptc}.boundary_list)
      continue
    end
    x = squeeze (msh{iptc}.boundary(iside).geo_map(1,:,:));
    y = squeeze (msh{iptc}.boundary(iside).geo_map(2,:,:));
    z = squeeze (msh{iptc}.boundary(iside).geo_map(3,:,:));
    gval = reshape (g (x, y, z, iref), ...
           msh{iptc}.boundary(iside).nqn, msh{iptc}.boundary(iside).nel);
    rhs_nmnn = ...
           op_f_v (sp{iptc}.boundary(iside), msh{iptc}.boundary(iside), gval);
    global_dofs = gnum{iptc}(sp{iptc}.boundary(iside).dofs);
    rhs(global_dofs) = rhs(global_dofs) + rhs_nmnn;
  end
end

% Apply Dirichlet boundary conditions
u = zeros (ndof, 1);
[u_drchlt, drchlt_dofs] = mp_mpi_sp_drchlt_l2_proj(mpi_comm, drchlt_prec, drchlt_solv, sp, msh, h, gnum, boundaries, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt;
rhs = MPI_Allreduce(rhs, 'mpi_sum', mpi_comm);

% Solve the linear system
preconditioner = mpi_prec( stiff_mat(int_dofs,int_dofs), mpi_comm );
u(int_dofs)    = mpi_solv( stiff_mat(int_dofs,int_dofs), rhs(int_dofs), preconditioner,  mpi_comm);

end

%!demo
%! ex_laplace_cube_mp

%!demo
%! ex_laplace_thick_L_mp