% EX_BSPLINE_STOKES_2D_STREAM: Solve a Stokes flow in stream-function formulation 
%                              problem on a two-dimensional domain.
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
geometry = geo_load (eye (4));

% Construct msh structure
breaks = {(linspace (0, 1, 11)), (linspace (0, 1, 11))};
knotsp = kntbrkdegreg (breaks, [3 3], [2 2]);
[qn, qw] = msh_set_quad_nodes (breaks, msh_gauss_nodes ([5 5]));
msh = msh_2d_tensor_product (breaks, qn, qw); 
msh = msh_push_forward_2d (msh, geometry, 'der2', true);

% Construct space structure
sp = sp_bspline_2d_phys (knotsp, [3 3], msh, 'hessian', true);

% Assemble the matrices
[x, y] = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)));
A = op_gradcurlu_gradcurlv_2d (sp, sp, msh, ones (size (x))); 

fx = @(x, y) (exp(x).*(2.*(-1 + y).*y.*(2 + (-5 + y).*y) - ...
              8.*x.*(-1 + y).*y.*(2 + (-5 + y).*y) + ...
              6.*x.^3.*(-4 + (-5 + y).*(-2 + y).*y.*(1 + y)) + ...
              x.^2.*(12 + y.*(-38 + y.*(19 + (-6 + y).*y))) + ...
              x.^4.*(12 + y.*(-38 + y.*(19 + (-6 + y).*y)))));

fy = @(x, y) (456 - 912.*y + ...
              exp(x).*(-456 - 2.*x.*(-23 + 5.*x).*(10 + (-3 + x).*x) + ...
              2.*(456 + x.*(-466 + x.*(253 + x.*(-82 + 7.*x)))).* ...
              y  + (-6 + x.*(6 + x.*(-11 + x.*(22 + 7.*x)))).* ...
              y.^2 + 2.*(6 + x.*(10 + x.*(-29 + (-6 + x).*x))).* ...
              y.^3 + (-6 + x.*(3 + x).*(-2 + x.*(7 + x))).*y.^4));

f = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));

b = op_f_curlv_2d (sp, msh, f (x, y));

% Apply Dirichlet boundary conditions
drchlt_dofs = unique ([sp.boundary(:).dofs]);
int_dofs = setdiff (1:sp.ndof, drchlt_dofs);


u = zeros (sp.ndof, 1);
u(int_dofs) = A(int_dofs, int_dofs) \ b(int_dofs);

sp_to_vtk_2d (u, sp, geometry, [61 61], 'stream_solution.vts', 'u')

pts = {linspace(0,1,31)', linspace(0,1,31)'};
[eu, F] = sp_eval_curl_2d (u, sp, geometry, pts);
msh_to_vtk_2d (F, squeeze (eu), 'vel.vts', 'vel')
