% EX_ARTICLE_SECTION_514: short example with non-homogeneous boundary conditions.
%
% Example to solve the problem
%
%    - div ( grad (u)) = 0                  in Omega
%                du/dn = -exp(x) * cos(y)   on Gamma_N
%                    u =  exp(x) * cos(y)   on Gamma_D
%
% with                Omega = (1 < x^2+y^2 < 4) & (x > 0) & (y > 0)
% and exact solution      u = exp(x) * cos(y)
%
% This solves the example of Section 5.1.4 in the article
%
% C. De Falco, A. Reali, R. Vazquez
% GeoPDEs: a research tool for IsoGeometric Analysis of PDEs
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

drchlt_sides = [1 2 4];
nmnn_sides   = 3;

geometry = geo_load ('ring_refined.mat');
knots = geometry.nurbs.knots;

[qn, qw] = msh_set_quad_nodes (knots, msh_gauss_nodes (geometry.nurbs.order));
msh = msh_2d_tensor_product (knots, qn, qw); 
msh = msh_push_forward_2d (msh, geometry);

space  = sp_nurbs_2d_phys (geometry.nurbs, msh);

[x, y] = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)));
mat = op_gradu_gradv (space, space, msh, ones (size (x))); 
rhs = zeros (space.ndof, 1);

for iside = nmnn_sides
  x = squeeze (msh.boundary(iside).geo_map(1,:,:));
  y = squeeze (msh.boundary(iside).geo_map(2,:,:));
  gval = - exp(x) .* cos(y);
  rhs_loc = op_f_v (space.boundary(iside), msh.boundary(iside), gval);
  rhs(space.boundary(iside).dofs) = rhs(space.boundary(iside).dofs) + rhs_loc;
end

drchlt_dofs = unique ([space.boundary(drchlt_sides).dofs]);
int_dofs = setdiff (1:space.ndof, drchlt_dofs);
M_drchlt = spalloc (space.ndof, space.ndof, space.ndof);
rhs_drchlt = zeros (space.ndof, 1);

for iside = drchlt_sides
  sp_bnd  = space.boundary(iside);
  msh_bnd = msh.boundary(iside);
  x = squeeze (msh_bnd.geo_map(1,:,:));
  y = squeeze (msh_bnd.geo_map(2,:,:));
  hval = exp(x) .* sin(y);
  M_side = op_u_v (sp_bnd, sp_bnd, msh_bnd, ones (size(x)));
  M_drchlt(sp_bnd.dofs, sp_bnd.dofs) = M_drchlt(sp_bnd.dofs, sp_bnd.dofs) + M_side;
  rhs_side = op_f_v (sp_bnd, msh_bnd, hval);
  rhs_drchlt(sp_bnd.dofs) = rhs_drchlt(sp_bnd.dofs) + rhs_side;
end

u = zeros (space.ndof, 1);
u(drchlt_dofs) = M_drchlt(drchlt_dofs, drchlt_dofs) \ rhs_drchlt(drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - mat(int_dofs, drchlt_dofs) * ...
    u(drchlt_dofs);

u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

sp_to_vtk_2d (u, space, geometry, [20 20], 'laplace_solution.vts', 'u')
err = sp_l2_error (space, msh, u, @(x,y) exp(x) .* sin(y))

%!demo
%! ex_article_section_514;
