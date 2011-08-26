% EX_ARTICLE_SECTION_513: short example with macro-element integration rule.
%
% Example to solve the problem
%
%    - div ( grad (u)) = (8-9*sqrt(x.^2+y.^2)).*sin(2*atan(y./x))./(x.^2+y.^2)  in Omega
%                    u = 0                                                      on Gamma
%
% with                Omega = (1 < x^2+y^2 < 4) & (x > 0) & (y > 0)
% and exact solution      u = (x.^2+y.^2-3*sqrt(x.^2+y.^2)+2).*sin(2.*atan(y./x))
%
% using a quadrature rule on macro-elements as proposed in
%
% 
% [1] Hughes, Reali, Sangalli. Comput. Methods Appl. Mech. Engrg. 
% Volume 199, Issues 5-8, 1 January 2010, Pages 301-313 
%
% This is the example of Section 5.1.3 in the article
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


nodes   = [5.168367524056075e-02, 2.149829914261059e-01, ...
           3.547033685486441e-01, 5.000000000000000e-01, ...
           6.452966314513557e-01, 7.850170085738940e-01, ...
           9.483163247594394e-01];
weights = [1.254676875668223e-01, 1.708286087294738e-01, ...
           1.218323586744639e-01, 1.637426900584793e-01, ...
           1.218323586744638e-01, 1.708286087294741e-01, ...
           1.254676875668223e-01];
rule    = [nodes; weights];

geometry = geo_load ({@ring_polar_map, @ring_polar_map_der});

[knots, breaks] = kntuniform ([10 10], [2 2], [1 1]);
breaks{1} ([2:3:end-2,3:3:end-1]) = [];
breaks{2} ([2:3:end-2,3:3:end-1]) = [];
[qn, qw]  = msh_set_quad_nodes (breaks, {rule, rule}, [0 1]);
msh = msh_2d (breaks, qn, qw, geometry);

space = sp_bspline_2d (knots, [2, 2], msh);

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
%! ex_article_section_513;
