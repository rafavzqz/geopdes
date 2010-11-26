% EX_BSPLINE_LAPLACE_3D: Solve a 3d Laplace problem with a B-spline discretization. 
%
% Example to solve the diffusion problem
%
% - div ( epsilon(x) grad (u)) = f    in Omega = F((0,1)^3)
%             epsilon(x) du/dn = g    on Gamma_N
%                            u = h    on Gamma_D
%
% In order to make it work, a script file with the data must also be executed.
% The available data files for this problem are:
%
% - test_cube
% - test_thick_ring
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
geometry = geo_load (geo_name);

[knots, zeta] = kntrefine (geometry.nurbs.knots, n_sub, degree, regularity);
  
% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_3d_tensor_product (zeta, qn, qw);
msh      = msh_push_forward_3d (msh, geometry);

% Construct space structure
sp       = sp_bspline_3d_phys (knots, degree, msh);

% Precompute the coefficients
x = squeeze (msh.geo_map(1,:,:));
y = squeeze (msh.geo_map(2,:,:));
z = squeeze (msh.geo_map(3,:,:));
  
epsilon = reshape (c_diff (x, y, z), msh.nqn, msh.nel);
fval    = reshape (f (x, y, z), msh.nqn, msh.nel) ;
 
% Assemble the matrices
stiff_mat = op_gradu_gradv (sp, sp, msh, epsilon);
rhs       = op_f_v (sp, msh, fval);

% Apply Neumann boundary conditions
for iside = nmnn_sides
  x = squeeze (msh.boundary(iside).geo_map(1,:,:));
  y = squeeze (msh.boundary(iside).geo_map(2,:,:));
  z = squeeze (msh.boundary(iside).geo_map(3,:,:));
  gval = reshape (g (x, y, z, iside), msh.boundary(iside).nqn, msh.boundary(iside).nel);

  rhs(sp.boundary(iside).dofs) = rhs(sp.boundary(iside).dofs) + ...
      op_f_v (sp.boundary(iside), msh.boundary(iside), gval);
end

% Apply Dirichlet boundary conditions
u = zeros (sp.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj(sp, msh, h, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:sp.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system
u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);

% Postprocessing
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk_3d (u, sp, geometry, vtk_pts, output_file, 'u');

if (exist ('uex', 'var'))
  error_l2 = sp_l2_error (sp, msh, u, uex)
  if (exist ('graduex', 'var'))
    error_h1 = sp_h1_error (sp, msh, u, uex, graduex)
  end
end

%!demo
%! test_thick_ring
%! ex_bspline_laplace_3d 
  
%!demo
%! test_cube
%! ex_bspline_laplace_3d 
