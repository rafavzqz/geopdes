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
geometry = geo_load ('annulus.mat');

% Construct msh structure
breaks = {(linspace (0, 1, 5)), (linspace (0, 1, 5))};
knotsp = kntbrkdegreg (breaks, [3 3], [2 2]);
[qn, qw] = msh_set_quad_nodes (breaks, msh_gauss_nodes ([5 5]));
msh = msh_2d_tensor_product (breaks, qn, qw); 
msh = msh_push_forward_2d (msh, geometry, 'der2', true);

% Construct space structure
sp = sp_bspline_2d_phys (knotsp, [3 3], msh, 'hessian', true);

% Assemble the matrices
[x, y] = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)));
A = op_gradcurlu_gradcurlv_2d (sp, sp, msh, ones (size (x))); 

fx = @(x, y) (4*(20*(98-57*x.^2).*y.^7+42*x.*(34*x.^2-75).*y.^6 ...
                 -3*(430*x.^4-1480*x.^2+1023).*y.^5 ...
                 +15*x.*(70*x.^4-310*x.^2+297).*y.^4 ...
                 +3*x.*(x.^4-5*x.^2+4).^2 ...
                 +2*(-290*x.^6+1500*x.^4-2079*x.^2+680).*y.^3 ...
                 +6*x.*(38*x.^6-255*x.^4+495*x.^2-260).*y.^2 ...
                 +(-75*x.^8+520*x.^6-1089*x.^4+720*x.^2-112).*y ...
                 +603*x.*y.^8-355*y.^9)-y./(x.^2+y.^2));

fy = @(x, y) (x./(x.^2+y.^2)+4*(5*x.^9-27*x.^8.*y+20*x.^7.*(15*y.^2-2) ...
                                +14*x.^6.*y.*(15-38*y.^2)+3*x.^5.*(290*y.^4-520*y.^2+33) ...
                                -15*x.^4.*y.*(70*y.^4-170*y.^2+33) ...
                                +x.^3.*(860*y.^6-3000*y.^4+2178*y.^2-80) ...
                                +18*x.^2.*y.*(-34*y.^6+155*y.^4-165*y.^2+20) ...
                                +x.*(285*y.^8-1480*y.^6+2079*y.^4-720*y.^2+16) ...
                                -67*y.^9+450*y.^7-891*y.^5+520*y.^3-48*y));
f = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));

b = op_f_curlv_2d (sp, msh, f (x, y));

% Apply Dirichlet boundary conditions
drchlt_dofs = unique ([sp.boundary(:).dofs]);
int_dofs = setdiff (1:sp.ndof, drchlt_dofs);


u = zeros (sp.ndof, 1);
u(int_dofs) = A(int_dofs, int_dofs) \ b(int_dofs);

sp_to_vtk_2d (u, sp, geometry, [20 20], 'stream_solution.vts', 'u')
