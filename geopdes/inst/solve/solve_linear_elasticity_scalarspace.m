% SOLVE_LINEAR_ELASTICITY_SCALARSPACE: Solve a linear elasticity problem on a NURBS domain.
%     It uses a scalar space of type sp_scalar, instead of sp_vector.
%
% The function solves the linear elasticity problem
%
%      - div (sigma(u)) = f    in Omega = F((0,1)^n)
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
%  [geometry, msh, space, u] = solve_linear_elasticity_scalarspace (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition (ONLY HOMOGENEOUS CONDITIONS)
%    - press_sides:  sides with pressure boundary condition (NOT IMPLEMENTED)
%    - symm_sides:   sides with symmetry boundary condition (NOT IMPLEMENTED)
%    - lambda_lame:  first Lame' parameter
%    - mu_lame:      second Lame' parameter
%    - f:            source term
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
%  msh:      mesh object that defines the quadrature rule (see msh_cartesian)
%  space:    space object that defines the discrete basis functions (see sp_scalar)
%  u:        the computed degrees of freedom
%
% See also EX_LIN_ELAST_HORSESHOE for an example.
%
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2011, 2015, 2022 Rafael Vazquez
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
              solve_linear_elasticity_scalarspace (problem_data, method_data)

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
msh      = msh_cartesian (geometry.nurbs.knots, qn, qw, geometry);

% Construct space structure
sp = sp_nurbs (nurbs, msh);

% Assemble the matrices
mat    = op_su_ev_tp (sp, sp, msh, lambda_lame, mu_lame); 
rhs    = op_f_v_tp_vector (sp, msh, f);

% Apply Neumann boundary conditions
for iside = nmnn_sides
% Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
  gside = @(varargin) g(varargin{:},iside);
  scalar_dofs = sp.boundary(iside).dofs;
  dofs = [];
  for icomp = 1:msh.rdim
    dofs = union (dofs, (icomp-1)*sp.ndof + scalar_dofs);
  end
  rhs_loc = op_f_v_tp_vector (sp.boundary(iside), msh.boundary(iside), gside);
  rhs(dofs) = rhs(dofs) + rhs_loc;
end

% % % Apply pressure conditions
% % for iside = press_sides
% %   msh_side = msh_eval_boundary_side (msh, iside);
% %   sp_side  = sp_eval_boundary_side (sp, msh_side);
% % 
% %   x = cell (msh_side.rdim, 1);
% %   for idim = 1:msh_side.rdim
% %     x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
% %   end
% %   pval = reshape (p (x{:}, iside), msh_side.nqn, msh_side.nel);
% % 
% %   rhs(sp_side.dofs) = rhs(sp_side.dofs) - op_pn_v (sp_side, msh_side, pval);
% % end
% % 
% % % Apply symmetry conditions
% % symm_dofs = [];
% % for iside = symm_sides
% %   if (~strcmpi (sp.transform, 'grad-preserving'))
% %     error ('The symmetry condition is only implemented for spaces with grad-preserving transform')
% %   end
% %   msh_side = msh_eval_boundary_side (msh, iside);
% %   for idim = 1:msh.rdim
% %     normal_comp(idim,:) = reshape (msh_side.normal(idim,:,:), 1, msh_side.nqn*msh_side.nel);
% %   end
% % 
% %   parallel_to_axes = false;
% %   for ind = 1:msh.rdim
% %     ind2 = setdiff (1:msh.rdim, ind);
% %     if (all (all (abs (normal_comp(ind2,:)) < 1e-10)))
% %       symm_dofs = union (symm_dofs, sp.boundary(iside).dofs(sp.boundary(iside).comp_dofs{ind}));
% %       parallel_to_axes = true;
% %       break
% %     end
% %   end
% %   if (~parallel_to_axes)
% %     error ('solve_linear_elasticity: We have only implemented the symmetry condition for boundaries parallel to the axes')
% %   end
% % 
% % end

% Apply Dirichlet boundary conditions
ndof = msh.rdim * sp.ndof;
u = zeros (ndof, 1);
% TO BE DONE. FOR NOW, ONLY HOMOGENEOUS CONDITIONS
drchlt_dofs = [];
for iside = 1:numel(drchlt_sides)
  side = drchlt_sides(iside);
  if (~exist('drchlt_components','var'))
    components = 1:msh.rdim;
  else
    components = drchlt_components{iside};
  end
  scalar_dofs = sp.boundary(side).dofs;
  for icomp = components
    drchlt_dofs = union (drchlt_dofs, (icomp-1)*sp.ndof + scalar_dofs);
  end
end

int_dofs = setdiff (1:ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - mat(int_dofs, drchlt_dofs)*u(drchlt_dofs);

% Solve the linear system
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

end
