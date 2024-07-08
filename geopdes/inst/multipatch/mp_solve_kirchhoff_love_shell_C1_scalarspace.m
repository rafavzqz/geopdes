% MP_SOLVE_KIRCHHOFF_LOVE_SHELL_C1_SCALARSPACE: Solve the Kirchhoff-Love shell model in an AS-G1 multipatch domain.
%  It uses an object space of type sp_multipatch_C1 with scalar values, instead of vector valued.
%
% USAGE:
%
%  [geometry, msh, space, u] = mp_solve_kirchhoff_love_shell_C1_scalarspace (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - drchlt_components: cell-array, the components that are set to zero for each drchlt_side
%    - E_coeff:      function handle for Young's modulus
%    - nu_coeff:     function handle for Poisson's ratio
%    - thickness:    scalar value, thickness of the shell
%    - f:            source term, distributed load
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree (>=3) of the spline functions.
%    - regularity: continuity (>=1, <=degree-2) of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%
% OUTPUT:
%
%  geometry: geometry structure (see mp_geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_multipatch)
%  space:    space object that defines the discrete basis functions (see sp_multipatch_C1)
%  u:        the computed degrees of freedom
%
% Copyright (C) 2017-2019 Pablo Antolin, Luca Coradello, Rafael Vazquez
% Copyright (C) 2022 Rafael Vazquez
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

function [geometry, msh, space, u] = ...
              mp_solve_kirchhoff_love_shell_C1_scalarspace (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

if (any(degree <= 1) || any(regularity == 0))
  error ('The degree must be at least two, and the regularity at least C^1')
end

% Construct geometry structure, and information for interfaces and boundaries
[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (geo_name);
npatch = numel (geometry);

msh = cell (1, npatch); 
sp = cell (1, npatch);
for iptc = 1:npatch

% Define the refined mesh, with tensor product structure
  [knots{iptc}, zeta{iptc}] = ...
         kntrefine (geometry(iptc).nurbs.knots, nsub-1, degree, regularity);

% Compute the quadrature rule
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (zeta{iptc}, rule);
  msh{iptc} = msh_cartesian (zeta{iptc}, qn, qw, geometry(iptc));

% Evaluate the discrete space basis functions in the quadrature points
  sp{iptc} = sp_bspline (knots{iptc}, degree, msh{iptc});
end

[edges, vertices] = vertices_struct (geometry, interfaces, boundaries, boundary_interfaces);
msh = msh_multipatch (msh, boundaries);
space = sp_multipatch_C1 (sp, msh, geometry, edges, vertices);
clear sp

% Compute and assemble the matrices
K = op_KL_shells_mp (space, space, msh, E_coeff, nu_coeff, thickness);
rhs = op_f_v_mp_vector (space, msh, f);
% Apply zero rotation with Nitsche method
if (exist ('rotation_sides', 'var') && ~isempty(rotation_sides))
  K = K + sp_nitsche_KL_rotation (space, msh, rotation_sides, E_coeff, nu_coeff, thickness, penalty_coeff);
end

% Apply boundary conditions
% if (~isfield (problem_data, 'uex'))
  [u_drchlt, drchlt_dofs, kernel_dofs] = sp_drchlt_C1_shells (space, msh, drchlt_sides, drchlt_components);
% else
%   [u_drchlt, drchlt_dofs, kernel_dofs] = sp_drchlt_C1_exact_shells (space, msh, drchlt_sides, problem_data.uex);
% end

ndof = msh.rdim * space.ndof;
u = zeros (ndof, 1);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:ndof, drchlt_dofs);
add_dofs = kernel_dofs.quasi_interior_dofs; %this will contain the "boundary" vertex dofs which have been removed from drchlt_dofs
add_dofs = [add_dofs, space.ndof+add_dofs, 2*space.ndof+add_dofs];

% We assemble the (pieces of the) stiffness matrix, the rhs (and its correction taking 
% into account the Dirichlet conditions), and the basis change matrix (we will need it 
% to go from the basis with kernel vectors obtained when examining the Dirichlet conditions 
% to the usual basis)
vertex_dofs = kernel_dofs.all_vertex_dofs;
vertex_dofs = [vertex_dofs, space.ndof+vertex_dofs, 2*space.ndof+vertex_dofs];
% B_change = speye (space.ndof); %basis change matrix
% B_change(kernel_dofs.all_vertex_dofs,kernel_dofs.quasi_interior_dofs) = kernel_dofs.B_change_local;

B_change_vector = blkdiag(kernel_dofs.B_change_local, kernel_dofs.B_change_local, kernel_dofs.B_change_local);
K(:,add_dofs) = K(:,vertex_dofs) * B_change_vector;
K(add_dofs,:) = B_change_vector.' * K(vertex_dofs,:);
rhs(add_dofs) = B_change_vector.' * rhs(vertex_dofs);

% rhs(int_dofs) = rhs(int_dofs) - K(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system
u(int_dofs) = K(int_dofs, int_dofs) \ rhs(int_dofs);

% Switching to the usual basis using the local matrix for the vertex dofs
u_old = u(setdiff(vertex_dofs, add_dofs)); % Coefficients of the vertex functions that were already in the old basis
u(vertex_dofs) = B_change_vector * u(add_dofs);
u(setdiff(vertex_dofs, add_dofs)) = u(setdiff(vertex_dofs, add_dofs)) + u_old;

end
