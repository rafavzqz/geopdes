% EX_ARTICLE_SECTION_511: short example that performs k-refinement.
%
% Example to solve the problem
%
%    - div ( grad (u)) = (8-9*sqrt(x.^2+y.^2)).*sin(2*atan(y./x))./(x.^2+y.^2)  in Omega
%                    u = 0                                                      on Gamma
%
% with                Omega = (1 < x^2+y^2 < 4) & (x > 0) & (y > 0)
% and exact solution      u = (x.^2+y.^2-3*sqrt(x.^2+y.^2)+2).*sin(2.*atan(y./x))
%
% including degree elevation and knot insertion for refinement.
%
% This solves the example of Section 5.1.1 in the article
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


geometry = geo_load ('ring_refined.mat');
nurbs    = geometry.nurbs;
degelev  = max ([3 3] - (nurbs.order-1), 0);
nurbs    = nrbdegelev (nurbs, degelev);
[rknots, zeta, nknots] = kntrefine (nurbs.knots, [1 1], nurbs.order-1, nurbs.order-2);
nurbs    = nrbkntins (nurbs, nknots);
geometry = geo_load (nurbs);
knots    = geometry.nurbs.knots;

[qn, qw] = msh_set_quad_nodes (knots, msh_gauss_nodes (geometry.nurbs.order));
msh = msh_2d_tensor_product (knots, qn, qw); 
msh = msh_push_forward_2d (msh, geometry);

space  = sp_nurbs_2d_phys (geometry.nurbs, msh);

[x, y] = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)));
mat = op_gradu_gradv (space, space, msh, ones (size (x))); 
rhs = op_f_v (space, msh, (8-9*sqrt(x.^2+y.^2)).*sin(2*atan(y./x))./(x.^2+y.^2));

drchlt_dofs = unique ([space.boundary(:).dofs]);
int_dofs = setdiff (1:space.ndof, drchlt_dofs);

u = zeros (space.ndof, 1);
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

sp_to_vtk_2d (u, space, geometry, [20 20], 'laplace_solution.vts', 'u')
err = sp_l2_error (space, msh, u, @(x,y)(x.^2+y.^2-3*sqrt(x.^2+y.^2)+2).*sin(2.*atan(y./x)))

%!demo
%! ex_article_section_511;
