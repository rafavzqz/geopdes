% EX_BSPLINE_LAPLACE_2D_EIG: Compute eigenvalues and eigenvectors for the Laplace operator.
%
% Example to solve the diffusion eigenvalue problem
%
%  - div ( epsilon(x) grad (u)) = lambda (mu(x) u)  in Omega = F((0,1)^2)
%              epsilon(x) du/dn = 0                 on Gamma_N
%                             u = 0                 on Gamma_D
%
% In order to make it work, a script file with the data must also be executed.
% The available data files for this problem are:
%
% - test_square_eig
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
geometry  = geo_load (geo_name);

[knots, zeta] = kntrefine (geometry.nurbs.knots, n_sub, degree, regularity);
  
% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_2d_tensor_product (zeta, qn, qw);
msh      = msh_push_forward_2d (msh, geometry);
  
% Construct space structure
sp       = sp_bspline_2d_phys (knots, degree, msh);
  
% Precompute the coefficients
x = squeeze (msh.geo_map(1,:,:));
y = squeeze (msh.geo_map(2,:,:));
  
epsilon = reshape (c_diff (x, y), msh.nqn, msh.nel);
mu      = reshape (c_react (x, y), msh.nqn, msh.nel);
 
% Assemble the matrices
stiff_mat = op_gradu_gradv (sp, sp, msh, epsilon);
mass_mat  = op_u_v (sp, sp, msh, mu);

% Apply homogeneous Dirichlet boundary conditions
drchlt_dofs = unique ([sp.boundary(drchlt_sides).dofs]);
int_dofs = setdiff (1:sp.ndof, drchlt_dofs);
  
% Solve the eigenvalue problem
eigf = zeros (sp.ndof, numel(int_dofs));
[eigf(int_dofs, :), eigv] = ...
    eig (full (stiff_mat(int_dofs, int_dofs)), full (mass_mat(int_dofs, int_dofs)));
eigv = diag (eigv);

if (all (imag (eigv) == 0))
  [eigv, perm] = sort (eigv);
elseif (any (abs (imag (eigv)) > 1e-9))
  error ('Complex eigenvalues appeared. I skip the postprocess.')
else
  warning ('Complex eigenvalues appeared, with small imaginary part. Only the real part is used for postprocessing')
  [eigv, perm] = sort ( real (eigv));
end

% Postprocessing
up = linspace(0,1,31)';
vp = linspace(0,1,31)';

% Plot of the 11th eigenfunction
subplot(1,2,1)
[eu, F] = sp_eval_2d (eigf(:, perm(11)), sp, geometry, {up vp});
surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), eu)
title ('Plot of the 11^{th} eigenfunction')

% Comparison with the exact eigenvalues
ndofs_1 = repmat ([1:sp.ndof_dir(1)-2], sp.ndof_dir(2)-2, 1);
ndofs_2 = repmat ([1:sp.ndof_dir(2)-2]', 1, sp.ndof_dir(1)-2);
exact = pi * sqrt (ndofs_1.^2 + ndofs_2.^2);
exact = sort (exact(:));
spectrum = sqrt (eigv) ./ exact;
subplot(1,2,2)
plot (linspace (0, 1, numel (spectrum)), spectrum, '*')
title ('Ratio of numerical to exact eigenvalues')

%!demo
%! test_square_eig
%! ex_bspline_laplace_2d_eig
