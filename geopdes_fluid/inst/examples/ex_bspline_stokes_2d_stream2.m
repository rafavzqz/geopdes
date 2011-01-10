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

% Construct geometry structure and set coefficients 

% problem on the parametric domain
geometry = geo_load (eye (4));
fx = @(x, y) (-(4 *(-1 + x).^2 + 16* (-1 + x).* x + 4* x.^2).* (-1 + y).^2.* y ...
              - (4* (-1 + x).^2 + 16* (-1 + x).* x + 4* x.^2).* (-1 + y).* y.^2 ...
              - 2* (-1 + x).^2.* x.^2.* (4* (-1 + y) + 2* y) - 2* (-1 + x).^2.* x.^2.* (2* (-1 + y) + 4* y)); 
fy = @(x, y) (-(-4* (-1 + x) - 8* x).* (-1 + y).^2.* y.^2 - (-8* (-1 + x) - 4* x).* (-1 +  y).^2.* y.^2 + 2 *...
              (-1 + x).^2.* x.* (2* (-1 + y).^2 + 8* (-1 + y).* y + 2* y.^2) + ...
              2* (-1 + x).* x.^2 .*(2* (-1 + y).^2 + 8* (-1 + y).* y + 2* y.^2)); 
      
vexx = @(x, y) (2*(-1 + x).^2.* x.^2 .*(-1 + y).^2 .*y + 2* (-1 + x).^2.* x.^2.* (-1 + y).* y.^2);
vexy = @(x, y) (-2* (-1 + x).^2.* x.* (-1 + y).^2.* y.^2 - 2* (-1 + x).* x.^2.* (-1 + y).^2.* y.^2);

% problem on a NURBS domain
% geometry = geo_load ('annulus.mat');
% fx = @(x, y) (4*(20*(98-57*x.^2).*y.^7+42*x.*(34*x.^2-75).*y.^6 ...
%                  -3*(430*x.^4-1480*x.^2+1023).*y.^5 ...
%                  +15*x.*(70*x.^4-310*x.^2+297).*y.^4 ...
%                  +3*x.*(x.^4-5*x.^2+4).^2 ...
%                  +2*(-290*x.^6+1500*x.^4-2079*x.^2+680).*y.^3 ...
%                  +6*x.*(38*x.^6-255*x.^4+495*x.^2-260).*y.^2 ...
%                  +(-75*x.^8+520*x.^6-1089*x.^4+720*x.^2-112).*y ...
%                  +603*x.*y.^8-355*y.^9)-y./(x.^2+y.^2));

% fy = @(x, y) (x./(x.^2+y.^2)+4*(5*x.^9-27*x.^8.*y+20*x.^7.*(15*y.^2-2) ...
%                                 +14*x.^6.*y.*(15-38*y.^2)+3*x.^5.*(290*y.^4-520*y.^2+33) ...
%                                 -15*x.^4.*y.*(70*y.^4-170*y.^2+33) ...
%                                 +x.^3.*(860*y.^6-3000*y.^4+2178*y.^2-80) ...
%                                 +18*x.^2.*y.*(-34*y.^6+155*y.^4-165*y.^2+20) ...
%                                 +x.*(285*y.^8-1480*y.^6+2079*y.^4-720*y.^2+16) ...
%                                 -67*y.^9+450*y.^7-891*y.^5+520*y.^3-48*y));

% vexx = @(x, y) (2*(x-y).*y.*(-4+x.^2+y.^2).*(-1+x.^2+y.^2).* ...
%                 (x.^5-8*y-2*x.^4.*y+20*y.^3-6*y.^5+2*x.^2.*y.*(5-4*y.^2)+x.^3.*(-5 + 6*y.^2)+ ...
%                  x.*(4+5*y.^2.*(-3+y.^2))));

% vexy = @(x, y) (-2*(x - y).*y.^2.*(-4 + x.^2 + y.^2).*(-1 + x.^2 + y.^2).* ...
%                 (4 + 5*x.^2.*(-3 + x.^2) + 2*x.*(5 - 2*x.^2).*y + ...
%                  (-5 + 6*x.^2).*y.^2 - 4*x.*y.^3 + y.^4));

f = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));

vex  = @(x, y) cat(1, ...
                   reshape (vexx (x,y), [1, size(x)]), ...
                   reshape (vexy (x,y), [1, size(x)]));


% Construct msh structure
breaks = {(linspace (0, 1, 31)), (linspace (0, 1, 31))};
knotsp = kntbrkdegreg (breaks, [4 4], [3 3]);
[qn, qw] = msh_set_quad_nodes (breaks, msh_gauss_nodes ([6 6]));
msh = msh_2d_tensor_product (breaks, qn, qw); 
msh = msh_push_forward_2d (msh, geometry, 'der2', true);

% Construct space structure
sp = sp_bspline_curl_2d_phys (knotsp, [4 4], msh);

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
figure (13)
title ('components of exact solution - computed solution')
subplot (1, 2, 1)
surf (X, Y, squeeze(euex(1,:,:))-squeeze(eu(1,:,:)))
title ('x component')
subplot (1, 2, 2)
surf (X, Y, squeeze(euex(2,:,:))-squeeze(eu(2,:,:)))
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
