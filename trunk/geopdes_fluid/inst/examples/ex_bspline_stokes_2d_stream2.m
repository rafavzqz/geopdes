% EX_BSPLINE_STOKES_2D_STREAM_2: Solve a Stokes flow in stream-function formulation 
%                                problem on a two-dimensional domain.
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
breaks = {(linspace (0, 1, 21)), (linspace (0, 1, 21))};
knotsp = kntbrkdegreg (breaks, [4 4], [3 3]);
[qn, qw] = msh_set_quad_nodes (breaks, msh_gauss_nodes ([6 6]));
msh = msh_2d_tensor_product (breaks, qn, qw); 
msh = msh_push_forward_2d (msh, geometry, 'der2', true);

% Construct space structure
sp = sp_bspline_curl_2d_phys (knotsp, [4 4], msh);

fx = @(x, y) (-(2*(y.^2.*(y-1)+y.*(y-1).^2)).*((x-1).^2+4*x.*(x-1)+x.^2) - ...
	      12*(x.^2.*(x-1).^2).*(2*y-1));
fy = @(x, y) (+(2*(x.^2.*(x-1)+x.*(x-1).^2)).*((y-1).^2+4*y.*(y-1)+y.^2) + ...
	      12*(y.^2.*(y-1).^2).*(2*x-1));;
f = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));

vexx = @(x, y) (2*(y.^2.*(y-1)+y.*(y-1).^2).*x.^2.*(x-1).^2);
vexy = @(x, y) (-2*(x.^2.*(x-1)+x.*(x-1).^2).*y.^2.*(y-1).^2);
vex  = @(x, y) cat(1, ...
                reshape (vexx (x,y), [1, size(x)]), ...
                reshape (vexy (x,y), [1, size(x)]));

% Assemble the matrices
[x, y] = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)));
A = op_gradu_gradv (sp, sp, msh, ones (size (x))); 
b = op_f_v (sp, msh, f (x, y));

% Apply Dirichlet boundary conditions
drchlt_dofs = unique ([sp.boundary(:).dofs]);
int_dofs = setdiff (1:sp.ndof, drchlt_dofs);

u = zeros (sp.ndof, 1);
u(int_dofs) = A(int_dofs, int_dofs) \ b(int_dofs);

l2_err = sp_l2_error (sp, msh, u, vex)

% Plot computed solution
vtk_pts = {linspace(0,1,31)', linspace(0,1,31)'};
[eu, F] = sp_eval_2d (u, sp, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
figure (1);
quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
axis equal
title('Computed solution')
figure (11)
title ('components of computed solution')
subplot (1, 2, 1)
surf (X, Y, squeeze(eu(1,:,:)))
title ('x component')
subplot (1, 2, 2)
surf (X, Y, squeeze(eu(2,:,:)))
title ('y component')

% Plot exact solution
figure (2);
euex = vex (X, Y);
quiver (X, Y, squeeze(euex(1,:,:)), squeeze(euex(2,:,:)))
axis equal
title('Exact solution')
figure (12)
title ('components of exact solution')
subplot (1, 2, 1)
surf (X, Y, squeeze(euex(1,:,:)))
title ('x component')
subplot (1, 2, 2)
surf (X, Y, squeeze(euex(2,:,:)))
title ('y component')

% compute L2 distance of the exact solution from the approximation space
M  = op_u_v (sp, sp, msh, ones (size (x))); 
b  = op_f_v (sp, msh, vex (x, y));
u2 = zeros (size (u));
u2(int_dofs) = M(int_dofs, int_dofs) \ b(int_dofs);

interperr_l2 = sp_l2_error (sp, msh, u2, vex)

% Plot project of exact solution on the approximation space
[eui, F] = sp_eval_2d (u2, sp, geometry, vtk_pts);
figure (3);
quiver (X, Y, squeeze(eui(1,:,:)), squeeze(eui(2,:,:)))
axis equal
title('projection')

% For FEM, res would be non-zero only at the boundary, why is it
% different here?
res = A * u2 - b;
vtk_pts = {linspace(0.3,1-.3,31)', linspace(0.3,1-.3,31)'};
[eui, F] = sp_eval_2d (res, sp, geometry, vtk_pts);
figure (4);
surf (X, Y, squeeze(eui(1,:,:)))
axis equal
figure (5);
surf (X, Y, squeeze(eui(2,:,:)))
axis equal
title('consistency error')

