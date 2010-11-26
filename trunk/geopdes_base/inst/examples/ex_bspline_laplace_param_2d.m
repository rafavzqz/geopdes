% EX_BSPLINE_LAPLACE_PARAM_2D: simple example in the parametric domain.
%
% Example to solve the particular problem
%
%    - div ( grad (u)) = 2*pi^2*sin(pi*x)*sin(pi*y)  in Omega = (0,1)^2
%                    u = 0                           on Gamma
%
% with exact solution    u = sin(pi*x)*sin(pi*y)
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


knots = {[0 0 0 0 1/4 1/2 3/4 1 1 1 1], [0 0 0 0 1/3 2/3 1 1 1 1]};
[qn, qw] = msh_set_quad_nodes (knots, msh_gauss_nodes ([3 3]));
msh = msh_2d_tensor_product (knots, qn, qw); 

sp  = sp_bspline_2d_param (knots, [3 3], msh);

[x, y] = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)));
mat = op_gradu_gradv (sp, sp, msh, ones (size (x))); 
rhs = op_f_v (sp, msh, 2*pi^2*sin(pi*x).*sin(pi*y));

drchlt_dofs = unique ([sp.boundary(:).dofs]);
int_dofs = setdiff (1:sp.ndof, drchlt_dofs);

u = zeros (sp.ndof, 1);
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

err = sp_l2_error (sp, msh, u, @(x,y) sin(pi*x).*sin(pi*y))

%!demo
%! ex_bspline_laplace_2d_15lines;
