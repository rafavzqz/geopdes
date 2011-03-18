% EX_BSPLINE_MIXED_FORM2_2D_EIG: Solve the 2d Maxwell eigenvalue problem with the second mixed formulation, and a B-spline discretization.
%
% Example to solve the problem
%    curl (1/mu(x) curl (u)) = lambda (epsilon(x) u)   in Omega = F((0,1)^2)
%          div (epsilon(x) u) = 0                       in Omega 
%       (1/mu(x) curl(u)) x n = 0                       on Gamma_N
%                       u x n = 0                       on Gamma_D
%
% with the variational mixed formulation
%
%  \int (epsilon(x) u v) + \int (curl(v) p) = 0,       \forall v \in H_0(curl),
%                          \int (curl(u) q) = -lambda \int (sqrt(1/mu(x)) p q), 
%                                                          \forall q \in L^2_0.
%
% USAGE:
%
%  [geometry, msh, space, sp_mul, eigv, eigf] = 
%                  solve_maxwell_eig_mixed2_2d (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_elec_perm:  electric permittivity (epsilon in the equation)
%    - c_magn_perm:  magnetic permeability (mu in the equation)
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - n_sub:      number of subdivisions for refinement.
%    - nquad:      number of points for Gaussian quadrature rule
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh structure (see msh_push_forward_2d)
%  space:    space structure (see sp_bspline_curl_transform_2d)
%  sp_mul:   space structure for the multiplier (see sp_bspline_3_forms)
%  eigv:     the computed eigenvalues
%  eigf:     degrees of freedom of the associated eigenfunctions
%
% See also EX_MAXWELL_EIG_MIXED2_SQUARE for an example
%
% Copyright (C) 2010, 2011 Rafael Vazquez
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

function [geometry, msh, sp, sp_mul, eigv, eigf] = ...
              solve_maxwell_eig_mixed2_2d (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} sprintf('= problem_data.(data_names{iopt});')]);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} sprintf('= method_data.(data_names{iopt});')]);
end

% Construct geometry structure 
geometry = geo_load (geo_name);

[knots, zeta] = kntrefine (geometry.nurbs.knots, n_sub, degree, regularity);
[knots_u1, knots_u2, knots_mul, degree1, degree2, degree_mul] = ...
                                                  knt_derham (knots, degree);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_2d_tensor_product (zeta, qn, qw);
msh      = msh_push_forward_2d (msh, geometry);

% Construct the space structures for the field and the Lagrange multiplier
sp = sp_bspline_curl_transform_2d (knots_u1, knots_u2, degree1, degree2, msh);
sp_mul = sp_bspline_3_forms (knots_mul, degree-1, msh);

% Precompute the coefficients
x = squeeze (msh.geo_map(1,:,:));
y = squeeze (msh.geo_map(2,:,:));

epsilon = reshape (c_elec_perm (x, y), msh.nqn, msh.nel);
mu      = reshape (c_magn_perm (x, y), msh.nqn, msh.nel);
coeff   = ones (msh.nqn, msh.nel);

% Assemble the matrices
mass_mat     = op_u_v (sp, sp, msh, epsilon);
mass_mul_mat = op_u_v (sp_mul, sp_mul, msh, sqrt(1./mu));
saddle_mat   = op_curlv_p (sp, sp_mul, msh, coeff);

% Apply homogeneous Dirichlet boundary conditions
drchlt_dofs  = unique ([sp.boundary(drchlt_sides).dofs]);
int_dofs     = setdiff (1:sp.ndof, drchlt_dofs);

% Solve the eigenvalue problem
mass_mat   = mass_mat (int_dofs, int_dofs);
saddle_mat = saddle_mat (:, int_dofs);

A = [mass_mat, saddle_mat.'; ...
     saddle_mat, sparse(sp_mul.ndof,sp_mul.ndof)];
M = [sparse(numel(int_dofs), numel(int_dofs) + sp_mul.ndof); ...
     sparse(sp_mul.ndof, numel(int_dofs)), -mass_mul_mat];

eigf = zeros (sp.ndof + sp_mul.ndof, numel(int_dofs) + sp_mul.ndof);
[eigf([int_dofs sp.ndof+[1:sp_mul.ndof]], :), eigv] = eig (full(A), full(M));
eigv = diag (eigv);

end