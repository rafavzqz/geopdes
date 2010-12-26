% EX_BSPLINE_MIXED_FORM1_2D_EIG: Solve the 2d Maxwell eigenvalue problem with a mixed formulation, and a B-spline discretization.
%
% Example to solve the problem
%
%    curl ( 1/mu(x) curl (u)) = lambda (epsilon(x) u)   in Omega = F((0,1)^2)
%          div (epsilon(x) u) = 0                       in Omega 
%       (1/mu(x) curl(u)) x n = 0                       on Gamma_N
%                       u x n = 0                       on Gamma_D
%
% with the variational mixed formulation
%
%    \int (1/mu(x) curl(u) curl(v)) + \int (epsilon(x) v grad(p)) 
%                = lambda \int (epsilon(x) u v),   \forall v \in H_0(curl),
%                                     \int (epsilon(x) u grad(q)) = 0,  
%                                                  \forall q \in H^1_0.
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

% Construct the space structures for the field and the Lagrange multiplier
sp = sp_bspline_curl_transform_2d (knots_u1, knots_u2, degree1, degree2, msh);
sp_mul = sp_bspline_2d_phys (knots, degree, msh);

% Precompute the coefficients
x = squeeze (msh.geo_map(1,:,:));
y = squeeze (msh.geo_map(2,:,:));

epsilon = reshape (c_elec_perm (x, y), msh.nqn, msh.nel);
mu      = reshape (c_magn_perm (x, y), msh.nqn, msh.nel);

% Assemble the matrices
stiff_mat  = op_curlu_curlv_2d (sp, sp, msh, 1./mu);
mass_mat   = op_u_v (sp, sp, msh, epsilon);
saddle_mat = op_v_gradp (sp, sp_mul, msh, epsilon);

% Apply homogeneous Dirichlet boundary conditions
drchlt_dofs     = unique ([sp.boundary(drchlt_sides).dofs]);
drchlt_dofs_mul = unique ([sp_mul.boundary(drchlt_sides).dofs]);
int_dofs     = setdiff (1:sp.ndof, drchlt_dofs);
int_dofs_mul = setdiff (1:sp_mul.ndof, drchlt_dofs_mul);

% Solve the eigenvalue problem
stiff_mat  = stiff_mat (int_dofs, int_dofs);
mass_mat   = mass_mat (int_dofs, int_dofs);
saddle_mat = saddle_mat (int_dofs_mul, int_dofs);

A = [stiff_mat, saddle_mat.'; ...
     saddle_mat, sparse(numel(int_dofs_mul),numel(int_dofs_mul))];
M = [mass_mat, sparse(numel(int_dofs), numel(int_dofs_mul)); ...
     sparse(numel(int_dofs_mul), numel(int_dofs)+numel(int_dofs_mul))];

eigf = zeros (sp.ndof + sp_mul.ndof, numel(int_dofs) + numel(int_dofs_mul));
[eigf([int_dofs sp.ndof+int_dofs_mul], :), eigv] = eig (full(A), full(M));
[eigv, perm] = sort (diag (eigv));

% Postprocessing
fprintf ('First computed eigenvalues: \n')
disp (eigv(1:5))

up = linspace(0,1,30)';
vp = linspace(0,1,30)';

figure
[eu, F] = sp_eval_2d (eigf(:,perm(8)), sp, geometry, {up vp});
quiver (squeeze(F(1,:,:)), squeeze(F(2,:,:)), ...
        squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
title ('8^{th} eigenfunction')

%!demo
%! test_maxwell_square_eig
%! ex_bspline_mixed_form1_2d_eig

%!demo
%! test_maxwell_Lshaped_eig
%! ex_bspline_mixed_form1_2d_eig

%!demo
%! test_maxwell_ring_1eighth_eig
%! ex_bspline_mixed_form1_2d_eig

%!demo
%! test_maxwell_curvedL_eig
%! ex_bspline_mixed_form1_2d_eig

