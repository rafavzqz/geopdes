% EX_BSPLINE_MAXWELL_2D: Solve a 2d Maxwell problem with a B-spline discretization.
%
% Example to solve the problem
%
%    curl ( epsilon(x) curl (u)) + u = f    in Omega = F((0,1)^2)
%           (epsilon(x) curl(u)) x n = g    on Gamma_N
%                              u x n = h    on Gamma_D
%
% In order to make it work, a script file with the data must also be executed.
% The available data files for this problem are:
%
% - test_maxwell_square
% - test_maxwell_Lshaped
% - test_maxwell_ring
% - test_maxwell_pacman
%
% Copyright (C) 2010 Rafael Vazquez
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

[knots, zeta] = kntrefine (geometry.nurbs.knots, n_sub, degree, regularity);
[knots_u1, knots_u2, degree1, degree2] = knt_derham (knots, degree);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_2d_tensor_product (zeta, qn, qw);
msh      = msh_push_forward_2d (msh, geometry);

% Construct space structure
sp = sp_bspline_curl_transform_2d (knots_u1, knots_u2, degree1, degree2, msh);

% Precompute the coefficients
x = squeeze (msh.geo_map(1,:,:));
y = squeeze (msh.geo_map(2,:,:));

epsilon = reshape (c_stiff (x, y), msh.nqn, msh.nel);
mu      = reshape (c_mass (x, y), msh.nqn, msh.nel);
fval    = reshape (f(x, y), 2, msh.nqn, msh.nel);

% Assemble the matrices
stiff_mat = op_curlu_curlv_2d (sp, sp, msh, epsilon);
mass_mat  = op_u_v (sp, sp, msh, mu);
rhs       = op_f_v (sp, msh, fval);

% Apply Neumann boundary conditions
for iside = nmnn_sides
  x = squeeze (msh.boundary(iside).geo_map(1,:,:));
  y = squeeze (msh.boundary(iside).geo_map(2,:,:));
  gval = reshape (g (x, y, iside), 2, msh.boundary(iside).nqn, msh.boundary(iside).nel);

  rhs(sp.boundary(iside).dofs) = rhs(sp.boundary(iside).dofs) + ...
    op_f_v (sp.boundary(iside), msh.boundary(iside), gval);
end

% Apply Dirichlet boundary conditions
u = zeros (sp.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj_uxn_2d (sp, msh, h, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:sp.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt ...
                              - mass_mat(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system
u(int_dofs) = (stiff_mat(int_dofs, int_dofs) + ...
                mass_mat(int_dofs, int_dofs)) \ rhs(int_dofs);

% Postprocessing
[eu, F] = sp_eval_2d (u, sp, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

if (exist ('uex', 'var'))
  eu2     = uex (X, Y);

  subplot(1,2,2)
  quiver (X, Y, squeeze(eu2(1,:,:)), squeeze(eu2(2,:,:)))
  axis equal
  title('Exact solution')
  subplot(1,2,1)

  error_l2 = sp_l2_error (sp, msh, u, uex)
  if (exist ('curluex', 'var'))
    error_hcurl = sp_hcurl_error (sp, msh, u, uex, curluex)
  end
end

quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
axis equal
title('Computed solution')

fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk_2d (u, sp, geometry, vtk_pts, output_file, 'u')

%!demo
%! test_maxwell_square
%! ex_bspline_maxwell_2d

%!demo
%! test_maxwell_Lshaped
%! ex_bspline_maxwell_2d

%!demo
%! test_maxwell_ring
%! ex_bspline_maxwell_2d

%!demo
%! test_maxwell_pacman
%! ex_bspline_maxwell_2d

