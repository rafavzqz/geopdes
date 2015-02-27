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


geometry = geo_load ('ring_refined.mat');
nurbs    = geometry.nurbs;
degelev  = max ([3 3] - (nurbs.order-1), 0);
nurbs    = nrbdegelev (nurbs, degelev);
[rknots, zeta, nknots] = kntrefine (nurbs.knots, [1 1], nurbs.order-1, nurbs.order-2);
nurbs    = nrbkntins (nurbs, nknots);
geometry = geo_load (nurbs);
knots    = geometry.nurbs.knots;

[qn, qw] = msh_set_quad_nodes (knots, msh_gauss_nodes (geometry.nurbs.order));
msh = msh_2d (knots, qn, qw, geometry);

space  = sp_nurbs_2d (geometry.nurbs, msh);

mat = op_gradu_gradv_tp (space, space, msh, @(x, y) ones (size (x))); 
rhs = op_f_v_tp (space, msh, @(x, y) (8-9*sqrt(x.^2+y.^2)).*sin(2*atan(y./x))./(x.^2+y.^2));

drchlt_dofs = [];
for iside = 1:4
  drchlt_dofs = union (drchlt_dofs, space.boundary(iside).dofs);
end
int_dofs = setdiff (1:space.ndof, drchlt_dofs);

u = zeros (space.ndof, 1);
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

sp_to_vtk (u, space, geometry, [20 20], 'laplace_solution.vts', 'u')
err = sp_l2_error (space, msh, u, @(x,y)(x.^2+y.^2-3*sqrt(x.^2+y.^2)+2).*sin(2.*atan(y./x)))

%!demo
%! ex_article_section_511;
