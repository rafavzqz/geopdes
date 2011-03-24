% SOLVE_LAPLACE_3D: Solve a 3d Laplace problem with a B-spline discretization (non-isoparametric approach). 
%
% The function solves the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega = F((0,1)^3)
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% USAGE:
%
%  [geometry, msh, space, u] = solve_laplace_3d (problem_data, method_data)
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
%    - n_sub:      number of subdivisions for refinement.
%    - nquad:      number of points for Gaussian quadrature rule
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh structure (see msh_push_forward_3d)
%  space:    space structure (see sp_bspline_3d_phys)
%  u:        the computed degrees of freedom
%
% See also EX_LAPLACE_CUBE for an example.
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011, Rafael Vazquez
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

function [geometry, msh, space, u] = ...
              solve_laplace_3d (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} sprintf('= problem_data.(data_names{iopt});')]);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} sprintf('= method_data.(data_names{iopt});')]);
end

% Construct geometry structure
geometry  = geo_load (geo_name);

[knots, zeta] = kntrefine (geometry.nurbs.knots, n_sub, degree, regularity);
  
% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_3d_tensor_product (zeta, qn, qw);
msh      = msh_push_forward_3d (msh, geometry);
  
% Construct space structure
space    = sp_bspline_3d_phys (knots, degree, msh);
  
% Precompute the coefficients
x = squeeze (msh.geo_map(1,:,:));
y = squeeze (msh.geo_map(2,:,:));
z = squeeze (msh.geo_map(3,:,:));
  
epsilon = reshape (c_diff (x, y, z), msh.nqn, msh.nel);
fval    = reshape (f (x, y, z), msh.nqn, msh.nel) ;
 
% Assemble the matrices
stiff_mat = op_gradu_gradv (space, space, msh, epsilon);
rhs       = op_f_v (space, msh, fval);

% Apply Neumann boundary conditions
for iside = nmnn_sides
  x = squeeze (msh.boundary(iside).geo_map(1,:,:));
  y = squeeze (msh.boundary(iside).geo_map(2,:,:));
  z = squeeze (msh.boundary(iside).geo_map(3,:,:));
  gval = reshape (g (x, y, z, iside), msh.boundary(iside).nqn, msh.boundary(iside).nel);

  rhs(space.boundary(iside).dofs) = rhs(space.boundary(iside).dofs) + ...
      op_f_v (space.boundary(iside), msh.boundary(iside), gval);
end

% Apply Dirichlet boundary conditions
u = zeros (space.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj(space, msh, h, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:space.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system
u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);

end
