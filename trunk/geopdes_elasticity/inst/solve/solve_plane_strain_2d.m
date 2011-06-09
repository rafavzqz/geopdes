% SOLVE_PLANE_STRAIN_2D: Solve a plane-strain problem on a two-dimensional domain.
%
% The function solves the plane strain problem
%
%      - div (sigma(u)) = f    in Omega = F((0,1)^2)
%      sigma(u) \cdot n = g    on Gamma_N
%                     u = h    on Gamma_D
%
% with   sigma(u) = mu*(grad(u) + grad(u)^t) + lambda*div(u)*I.
%
%   u:          displacement vector
%   sigma:      Cauchy stress tensor
%   lambda, mu: Lame' parameters
%   I:          identity tensor
%
% USAGE:
%
%  [geometry, msh, space, u] = solve_plane_strain_2d (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - press_sides:  sides with pressure boundary condition (may be empty)
%    - symm_sides:   sides with symmetry boundary condition (may be empty)
%    - lam:          first Lame' parameter
%    - mu:           second Lame' parameter
%    - f:            source term
%    - h:            function for Dirichlet boundary condition
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - p:            pressure function (if press_sides is not empty)
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh structure (see msh_push_forward_3d)
%  space:    space structure (see sp_bspline_3d_phys)
%  u:        the computed degrees of freedom
%
% See also EX_PLANE_STRAIN_RING for an example.
%
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2011 Rafael Vazquez
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

function [geometry, msh, sp, u] = ...
              solve_plane_strain_2d (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Construct geometry structure
geometry = geo_load (geo_name);
nurbs    = geometry.nurbs;
degelev  = max (degree - (nurbs.order-1), 0);
nurbs    = nrbdegelev (nurbs, degelev);
[rknots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);

nurbs = nrbkntins (nurbs, nknots);
geometry = geo_load (nurbs);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (geometry.nurbs.knots, rule);
msh      = msh_2d_tensor_product (geometry.nurbs.knots, qn, qw);
msh      = msh_push_forward_2d (msh, geometry);

% Construct space structure
sp_scalar = sp_nurbs_2d_phys (nurbs, msh);
sp = sp_scalar_to_vector_2d (sp_scalar, sp_scalar, msh, 'divergence', true);

% Assemble the matrices
[x, y]    = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)));
coeff_lam = reshape (lam (x, y), msh.nqn, msh.nel);
coeff_mu  = reshape (mu (x, y), msh.nqn, msh.nel);
fval      = reshape (f (x, y), 2, msh.nqn, msh.nel);

mat       = op_su_ev (sp, sp, msh, coeff_lam, coeff_mu); 
rhs       = op_f_v (sp, msh, fval);

% Apply Neumann boundary conditions
for iside = nmnn_sides
  x = squeeze (msh.boundary(iside).geo_map(1,:,:));
  y = squeeze (msh.boundary(iside).geo_map(2,:,:));
  rhs(sp.boundary(iside).dofs) = rhs(sp.boundary(iside).dofs) + op_f_v (sp.boundary(iside), msh.boundary(iside), g (x, y, iside));
end

% Apply pressure conditions
for iside = press_sides
  x = squeeze (msh.boundary(iside).geo_map(1,:,:));
  y = squeeze (msh.boundary(iside).geo_map(2,:,:));
  rhs(sp.boundary(iside).dofs) = rhs(sp.boundary(iside).dofs) - op_pn_v (sp.boundary(iside), msh.boundary(iside), p (x, y, iside));
end

% Apply symmetry conditions
u = zeros (sp.ndof, 1);
symm_dofs = [];
for iside = symm_sides
  if (all (abs (msh.boundary(iside).normal(1,:,:)) < 1e-10))
    symm_dofs = unique ([symm_dofs sp.boundary(iside).comp_dofs{2}]);
  elseif (all (abs (msh.boundary(iside).normal(2,:,:)) < 1e-10))
    symm_dofs = unique ([symm_dofs sp.boundary(iside).comp_dofs{1}]);
  else
    error ('ex_nurbs_plane_strain_2d: We have only implemented the symmetry condition for boundaries parallel to the axes')
  end
end

% Apply Dirichlet boundary conditions
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (sp, msh, h, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:sp.ndof, [symm_dofs drchlt_dofs]);
rhs(int_dofs) = rhs(int_dofs) - mat(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

end

%!demo
%! ex_plane_strain_square

%!demo
%! ex_plane_strain_square_mixed_bc

%!demo
%! ex_plane_strain_plate

%!demo
%! ex_plane_strain_ring
