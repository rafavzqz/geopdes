% SOLVE_LINEAR_ELASTICITY_3D: Solve a linear elasticity problem on a 3-dimensional domain.
%
% The function solves the linear elasticity problem
%
%      - div (sigma(u)) = f    in Omega = F((0,1)^3)
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
%  [geometry, msh, space, u] = solve_linear_elasticity_3d (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - press_sides:  sides with pressure boundary condition (may be empty)
%    - symm_sides:   sides with symmetry boundary condition (may be empty)
%    - lambda_lame:  first Lame' parameter
%    - mu_lame:      second Lame' parameter
%    - f:            source term
%    - h:            function for Dirichlet boundary condition
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
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
%  msh:      mesh object that defines the quadrature rule (see msh_3d)
%  space:    space object that defines the discrete basis functions (see sp_vector_3d)
%  u:        the computed degrees of freedom
%
% See also EX_LIN_ELAST_HORSESHOE for an example.
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
              solve_linear_elasticity_3d (problem_data, method_data)

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
degelev  = max (degree - (geometry.nurbs.order-1), 0);
nurbs    = nrbdegelev (geometry.nurbs, degelev);
[rknots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);

nurbs    = nrbkntins (nurbs, nknots);
geometry = geo_load (nurbs);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (geometry.nurbs.knots, rule);
msh      = msh_3d (geometry.nurbs.knots, qn, qw, geometry);

% Construct space structure
sp_scalar = sp_nurbs_3d (nurbs, msh);
sp = sp_vector_3d (sp_scalar, sp_scalar, sp_scalar, msh);
clear sp_scalar

% Assemble the matrices
mat    = op_su_ev_tp (sp, sp, msh, lambda_lame, mu_lame); 
rhs    = op_f_v_tp (sp, msh, f);

% Apply Neumann boundary conditions
for iside = nmnn_sides
  msh_side = msh_eval_boundary_side (msh, iside);
  sp_side  = sp_eval_boundary_side (sp, msh_side);

  x = squeeze (msh_side.geo_map(1,:,:));
  y = squeeze (msh_side.geo_map(2,:,:));
  z = squeeze (msh_side.geo_map(3,:,:));
  gval = reshape (g (x, y, z, iside), 3, msh_side.nqn, msh_side.nel);
  rhs(sp.boundary(iside).dofs) = rhs(sp.boundary(iside).dofs) + op_f_v (sp_side, msh_side, gval);
end

% Apply pressure conditions
for iside = press_sides
  msh_side = msh_eval_boundary_side (msh, iside);
  sp_side  = sp_eval_boundary_side (sp, msh_side);

  x = squeeze (msh_side.geo_map(1,:,:));
  y = squeeze (msh_side.geo_map(2,:,:));
  z = squeeze (msh_side.geo_map(2,:,:));
  pval = reshape (p (x, y, z, iside), msh_side.nqn, msh_side.nel);

  rhs(sp_side.dofs) = rhs(sp_side.dofs) - op_pn_v (sp_side, msh_side, pval);
end

% Apply symmetry conditions
u = zeros (sp.ndof, 1);
symm_dofs = [];
for iside = symm_sides
  msh_side = msh_eval_boundary_side (msh, iside);
  normal_comp1 = reshape (msh_side.normal(1,:,:), 1, msh_side.nqn*msh_side.nel);
  normal_comp2 = reshape (msh_side.normal(2,:,:), 1, msh_side.nqn*msh_side.nel);
  normal_comp3 = reshape (msh_side.normal(3,:,:), 1, msh_side.nqn*msh_side.nel);
  if ((all (abs (normal_comp2) < 1e-10)) && ...
      (all (abs (normal_comp3) < 1e-10)))
    symm_dofs = union (symm_dofs, sp.boundary(iside).comp_dofs{1});
  elseif ((all (abs (normal_comp1) < 1e-10)) && ...
          (all (abs (normal_comp3) < 1e-10)))
    symm_dofs = union (symm_dofs, sp.boundary(iside).comp_dofs{2});
  elseif ((all (abs (normal_comp1) < 1e-10)) && ...
          (all (abs (normal_comp2) < 1e-10)))
    symm_dofs = union (symm_dofs, sp.boundary(iside).comp_dofs{3});
  else
    error ('ex_nurbs_linear_elasticity_3d: We have only implemented the symmetry condition for boundaries parallel to the axes')
  end
end

% Apply Dirichlet boundary conditions
u = zeros (sp.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (sp, msh, h, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:sp.ndof, [drchlt_dofs, symm_dofs]);
rhs(int_dofs) = rhs(int_dofs) - mat (int_dofs, drchlt_dofs) * u_drchlt;

% Solve the linear system
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

end
