% EX_NURBS_LINEAR_ELASTICITY_3D: Solve a linear elasticity problem on a 3-dimensional domain.
%
% Example to solve the linear elasticity problem
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
% In order to make it work, a script file with the data must also be executed.
% The available data files for this problem are:
%
% - test_linear_elasticity_horseshoe
%
% Copyright (C) 2010 Carlo de Falco
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
geometry = geo_load (geo_name);
degelev  = max (degree - (geometry.nurbs.order-1), 0);
nurbs    = nrbdegelev (geometry.nurbs, degelev);
[rknots, zeta, nknots] = kntrefine (nurbs.knots, n_sub, nurbs.order-1, regularity);

nurbs    = nrbkntins (nurbs, nknots);
geometry = geo_load (nurbs);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (geometry.nurbs.knots, rule);
msh      = msh_3d_tensor_product (geometry.nurbs.knots, qn, qw);
msh      = msh_push_forward_3d (msh, geometry);

% Construct space structure
sp_scalar = sp_nurbs_3d_phys (nurbs, msh);
sp = sp_scalar_to_vector_3d (sp_scalar, sp_scalar, sp_scalar, msh, 'divergence', true);

% Assemble the matrices
[x, y, z] = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)), squeeze (msh.geo_map(3,:,:)));
mat    = op_su_ev (sp, sp, msh, lam (x, y, z), mu (x, y, z)); 
rhs    = op_f_v (sp, msh, f (x, y, z));

% Apply Neumann boundary conditions
for iside = nmnn_sides
  x = squeeze (msh.boundary(iside).geo_map(1,:,:));
  y = squeeze (msh.boundary(iside).geo_map(2,:,:));
  z = squeeze (msh.boundary(iside).geo_map(3,:,:));
  rhs(sp.boundary(iside).dofs) = rhs(sp.boundary(iside).dofs) + op_f_v (sp.boundary(iside), msh.boundary(iside), g (x, y, z, iside));
end

% Apply Dirichlet boundary conditions
u = zeros (sp.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (sp, msh, h, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:sp.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - mat (int_dofs, drchlt_dofs) * u_drchlt;

% Solve the linear system
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

% Postprocessing
if (exist ('uex', 'var'))
  l2_error = sp_l2_error (sp, msh, u, uex)
end

fprintf ('results being saved in: %s_displacement\n', output_file)
sp_to_vtk_3d (u, sp, geometry, vtk_pts, sprintf ('%s_displacement', output_file), 'displacement')

%!demo
%! test_linear_elasticity_horseshoe
%! ex_nurbs_linear_elasticity_3d
