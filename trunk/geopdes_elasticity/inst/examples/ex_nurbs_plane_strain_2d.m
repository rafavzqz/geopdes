% EX_NURBS_PLANE_STRAIN_2D: Solve a plane-strain problem on a two-dimensional domain.
%
% Example to solve the plane strain problem
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
% In order to make it work, a script file with the data must also be executed.
% The available data files for this problem are:
%
% - test_plane_strain_square
% - test_plane_strain_square_mixed_bc
% - test_plane_strain_plate
% - test_plane_strain_ring
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
nurbs    = geometry.nurbs;
degelev  = max (degree - (nurbs.order-1), 0);
nurbs    = nrbdegelev (nurbs, degelev);
[rknots, zeta, nknots] = kntrefine (nurbs.knots, n_sub, nurbs.order-1, regularity);

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
[x, y] = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)));
mat    = op_su_ev (sp, sp, msh, lam (x, y), mu (x, y)); 
rhs    = op_f_v (sp, msh, f (x, y));

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
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj(sp, msh, h, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:sp.ndof, [symm_dofs drchlt_dofs]);
rhs(int_dofs) = rhs(int_dofs) - mat(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

% Postprocessing
[eu, F] = sp_eval_2d (u, sp, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

figure()
if (exist ('uex', 'var'))
  l2_error = sp_l2_error (sp, msh, u, uex)
  
  subplot(1,2,2)
  eu2 = uex (X, Y);
  quiver (X, Y, squeeze(eu2(1,:,:)), squeeze(eu2(2,:,:)))
  axis equal
  title('Exact solution')
  
  subplot(1,2,1)
end

quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
axis equal
title('Computed solution')

fprintf ('results being saved in: %s_displacement\n', output_file)
sp_to_vtk_2d (u, sp, geometry, vtk_pts, sprintf ('%s_displacement.vts', output_file), 'displacement')

%!demo
%! test_plane_strain_square
%! ex_nurbs_plane_strain_2d

%!demo
%! test_plane_strain_square_mixed_bc
%! ex_nurbs_plane_strain_2d

%!demo
%! test_plane_strain_plate
%! ex_nurbs_plane_strain_2d

%!demo
%! test_plane_strain_ring
%! ex_nurbs_plane_strain_2d

