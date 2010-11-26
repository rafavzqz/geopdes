% EX_BSPLINE_MAXWELL_2D_EIG: Solve the 2d Maxwell eigenvalue problem with a B-spline discretization.
%
% Example to solve the problem
%
%    curl ( epsilon(x) curl (u)) = lambda (mu(x) u)   in Omega = F((0,1)^2)
%       (epsilon(x) curl(u)) x n = 0                  on Gamma_N
%                          u x n = 0                  on Gamma_D
%
% In order to make it work, a script file with the data must also be executed.
% The available data files for this problem are:
%
% - test_maxwell_square_eig
% - test_maxwell_Lshaped_eig
% - test_maxwell_ring_1eighth_eig
% - test_maxwell_curvedL_eig
%
% Copyright (C) 2010 Rafael Vazquez
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
geometry = geo_load (geo_name);

[knots, zeta] = kntrefine (geometry.nurbs.knots, n_sub, degree, regularity);
[knots_u1, knots_u2, degree1, degree2] = knt_derham (knots, degree);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_2d_tensor_product (zeta, qn, qw);
msh      = msh_push_forward_2d (msh, geometry);

% Construct space structure
sp = sp_bspline_curl_transform_2d (knots_u1, knots_u2, degree1, degree2, msh);

% Precompute the coefficients
x = squeeze (msh.geo_map(1,:,:));
y = squeeze (msh.geo_map(2,:,:));

epsilon = reshape (c_elec_perm (x, y), msh.nqn, msh.nel);
mu      = reshape (c_magn_perm (x, y), msh.nqn, msh.nel);

% Assemble the matrices
stiff_mat = op_curlu_curlv_2d (sp, sp, msh, 1./mu);
mass_mat  = op_u_v (sp, sp, msh, epsilon);

% Apply homogeneous Dirichlet boundary conditions
drchlt_dofs = unique ([sp.boundary(drchlt_sides).dofs]);
int_dofs = setdiff (1:sp.ndof, drchlt_dofs);

% Solve the eigenvalue problem
eigf = zeros (sp.ndof, numel(int_dofs));
[eigf(int_dofs, :), eigv] = ...
    eig (full (stiff_mat(int_dofs, int_dofs)), full (mass_mat(int_dofs, int_dofs)));
[eigv, perm] = sort (diag (eigv));


% Postprocessing
nzeros = numel (find (eigv < 1e-10));

fprintf ('Number of zero eigenvalues: %i \n', nzeros)
fprintf ('First nonzero eigenvalues: \n')
disp (eigv(nzeros+1:nzeros+5))

up = linspace(0,1,30)';
vp = linspace(0,1,30)';

figure
[eu, F] = sp_eval_2d (eigf(:,perm(nzeros+8)), sp, geometry, {up vp});
quiver (squeeze(F(1,:,:)), squeeze(F(2,:,:)), ...
        squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
title ('8^{th} eigenfunction')

%!demo
%! test_maxwell_square_eig
%! ex_bspline_maxwell_2d_eig

%!demo
%! test_maxwell_Lshaped_eig
%! ex_bspline_maxwell_2d_eig

%!demo
%! test_maxwell_ring_1eighth_eig
%! ex_bspline_maxwell_2d_eig

%!demo
%! test_maxwell_curvedL_eig
%! ex_bspline_maxwell_2d_eig

