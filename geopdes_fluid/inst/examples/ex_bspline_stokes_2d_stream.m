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
knotsp = kntbrkdegreg (breaks, [4 4], [3 3]);
[qn, qw] = msh_set_quad_nodes (breaks, msh_gauss_nodes ([5 5]));
msh = msh_2d_tensor_product (breaks, qn, qw); 
msh = msh_push_forward_2d (msh, geometry, 'der2', true);

% Construct space structure
sp = sp_bspline_2d_phys (knotsp, [4 4], msh, 'hessian', true);

% Assemble the matrices
[x, y] = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)));
A = op_gradcurlu_gradcurlv_2d (sp, sp, msh, ones (size (x))); 
%A2 = op_gradcurlu_gradcurlv_2d_bak (sp, sp, msh, ones (size (x))); 

fx = @(x, y) (-4.*(y.*(y-1).*2+y.^2.*(y-1)).*((x-1).^2+2.*x.*(x-1)+x.^2));
fy = @(x, y) (4.*(x.*(x-1).*2+x.^2.*(x-1)).*((y-1).^2+2.*y.*(y-1)+y.^2));

uex = @(x, y) (x.^2.*(x-1).^2.*y.^2.*(y-1).^2);

f = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));

b = op_f_curlv_2d (sp, msh, f (x, y));
%b2 = op_f_curlv_2d_bak (sp, msh, f (x, y));

M = op_u_v (sp, sp, msh, ones (size (x))); 
E = sum (M);

% Apply Dirichlet boundary conditions
drchlt_dofs = unique ([sp.boundary(:).dofs]);
int_dofs = setdiff (1:sp.ndof, drchlt_dofs);

u = zeros (sp.ndof, 1);
u(int_dofs) = A(int_dofs, int_dofs) \ b(int_dofs);


toterr = sp_l2_error (sp, msh, u, uex)

%sp_to_vtk_2d (u, sp, geometry, [61 61], 'stream_solution.vts', 'u')

pts = {linspace(0,1,31)', linspace(0,1,31)'};

[eu, F] = sp_eval_2d (u, sp, geometry, pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

figure(1);%subplot (1, 2, 1)
surf (X, Y, eu)
figure(2);%subplot (1, 2, 2)
surf (X, Y, uex (X, Y))

%[eu, F] = sp_eval_curl_2d (u, sp, geometry, pts);
%msh_to_vtk_2d (F, squeeze (eu), 'vel.vts', 'vel')

AA = op_u_v (sp, sp, msh, ones (size (x))); 
bb = op_f_v (sp, msh, uex (x, y)); 
unumex = AA \ bb;

interperr = sp_l2_error (sp, msh, unumex, uex)

[eu, F] = sp_eval_2d (unumex, sp, geometry, pts);

figure(3);%subplot (1, 2, 2)
surf (X, Y, eu)


res = A * unumex - b;
[eu, F] = sp_eval_2d (res, sp, geometry, pts);

figure(4);%subplot (1, 2, 2)
surf (X, Y, eu)

w = zeros (size (u));
w(int_dofs) = A(int_dofs, int_dofs) \ (b(int_dofs) - A(int_dofs, drchlt_dofs) * res(drchlt_dofs));
[eu, F] = sp_eval_2d (w, sp, geometry, pts);

figure(5);%subplot (1, 2, 2)
surf (X, Y, eu)
