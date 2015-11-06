% SOLVE_MAXWELL_EIG_MIXED2_2D: Solve the 2d Maxwell eigenvalue problem with the second mixed formulation, and a B-spline discretization.
%
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
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_2d)
%  space:    space object that defines the discrete functions (see sp_vector_2d_curl_transform)
%  sp_mul:   space object for the multiplier (see sp_bspline_2d_3forms)
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

function [geometry, msh, space, sp_mul, eigv, eigf] = ...
              solve_maxwell_eig_mixed2_2d (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Construct geometry structure 
geometry = geo_load (geo_name);

[knots, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);
[knots_u1, knots_u2, knots_mul, degree1, degree2, degree_mul] = ...
                                                  knt_derham (knots, degree);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_2d (zeta, qn, qw, geometry);

% Construct the space structures for the field and the Lagrange multiplier
sp_u1 = sp_bspline_2d (knots_u1, degree1, msh);
sp_u2 = sp_bspline_2d (knots_u2, degree2, msh);
space = sp_vector_2d_curl_transform (sp_u1, sp_u2, msh);
clear sp_u1 sp_u2
sp_mul = sp_bspline_2d_3forms (knots_mul, degree-1, msh);

% Assemble the matrices
sqrt_invmu = @(x, y) sqrt (1./c_magn_perm (x, y));
mass_mat  = op_u_v_tp (space, space, msh, c_elec_perm);
mass_mul_mat = op_u_v_tp (sp_mul, sp_mul, msh, sqrt_invmu);
saddle_mat   = op_curlv_p_tp (space, sp_mul, msh, @(x, y) ones (size (x)));

% Apply homogeneous Dirichlet boundary conditions
drchlt_dofs = [];
for iside = 1:numel (drchlt_sides)
  drchlt_dofs = union (drchlt_dofs, space.boundary(drchlt_sides(iside)).dofs);
end
int_dofs = setdiff (1:space.ndof, drchlt_dofs);

% Solve the eigenvalue problem
mass_mat   = mass_mat (int_dofs, int_dofs);
saddle_mat = saddle_mat (:, int_dofs);

A = [mass_mat, saddle_mat.'; ...
     saddle_mat, sparse(sp_mul.ndof,sp_mul.ndof)];
M = [sparse(numel(int_dofs), numel(int_dofs) + sp_mul.ndof); ...
     sparse(sp_mul.ndof, numel(int_dofs)), -mass_mul_mat];

eigf = zeros (space.ndof + sp_mul.ndof, numel(int_dofs) + sp_mul.ndof);
[eigf([int_dofs space.ndof+[1:sp_mul.ndof]], :), eigv] = eig (full(A), full(M));
eigv = diag (eigv);

end

%!demo
%! ex_maxwell_eig_mixed2_square

%!demo
%! ex_maxwell_eig_mixed2_Lshaped
