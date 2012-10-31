% EX_TSPLINES_15LINES: minimalistic T-splines tutorial in 15 lines.
%
% Example to solve the problem
%
%    - div ( grad (u)) = (8-9*sqrt(x.^2+y.^2)).*sin(2*atan(y./x))./(x.^2+y.^2)  in Omega
%                    u = 0                                                      on Gamma
%
% with                Omega = (1 < x^2+y^2 < 4) & (x > 0) & (y > 0)
% and exact solution      u = (x.^2+y.^2-3*sqrt(x.^2+y.^2)+2).*sin(2.*atan(y./x))
%
% with a T-spline discretization.
%
% This solves the example of Section 4 in the article
%
% C. De Falco, A. Reali, R. Vazquez
% GeoPDEs: a research tool for IsoGeometric Analysis of PDEs
%
% using a mesh with T-junctions.
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2012 Rafael Vazquez
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

tspline = read_bezier_extraction ('arch_tsplines.iga');
[msh, space] = tspline_mesh_space (tspline);

[x, y] = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)));
mat = op_gradu_gradv (space, space, msh, ones (size (x)));
rhs = op_f_v (space, msh, (8-9*sqrt(x.^2+y.^2)).*sin(2*atan(y./x))./(x.^2+y.^2));

drchlt_dofs = [];
for ibnd = 1:numel (tspline.boundary)
  [msh_bnd, sp_bnd] = tspline_boundary (tspline, ibnd);
  drchlt_dofs = union (drchlt_dofs, sp_bnd.connectivity(:)');
end
int_dofs = setdiff (1:space.ndof, drchlt_dofs);

u = zeros (space.ndof, 1);
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

tspline_to_vtk (u, tspline, [10 10], 'laplace_tspline_solution.vtu', 'u')
err = tspline_l2_error (space, msh, u, @(x,y)(x.^2+y.^2-3*sqrt(x.^2+y.^2)+2).*sin(2.*atan(y./x)))

%!demo
%! ex_tsplines_15lines;
