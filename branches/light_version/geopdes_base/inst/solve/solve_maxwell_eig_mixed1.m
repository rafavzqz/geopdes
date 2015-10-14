% SOLVE_MAXWELL_EIG_MIXED1: Solve the Maxwell eigenvalue problem with a mixed formulation, and a B-spline discretization.
%
% Example to solve the problem
%
%    curl (1/mu(x) curl (u)) = lambda (epsilon(x) u)   in Omega = F((0,1)^n)
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
% USAGE:
%
%  [geometry, msh, space, sp_mul, eigv, eigf] = 
%                  solve_maxwell_eig_mixed1 (problem_data, method_data)
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
%  msh:      mesh object that defines the quadrature rule (see msh_cartesian)
%  space:    space object that defines the discrete functions (see sp_vector_curl_transform)
%  sp_mul:   space object for the multiplier (see sp_bspline)
%  eigv:     the computed eigenvalues
%  eigf:     degrees of freedom of the associated eigenfunctions
%
% See also EX_MAXWELL_EIG_MIXED1_SQUARE for an example
%
% Copyright (C) 2010, 2011, 2015 Rafael Vazquez
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
              solve_maxwell_eig_mixed1 (problem_data, method_data)

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
[knots_hcurl, degree_hcurl] = knt_derham (knots, degree, 'Hcurl');

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);

% Construct the space structures for the field and the Lagrange multiplier
scalar_spaces = cell (msh.ndim, 1);
for idim = 1:msh.ndim
  scalar_spaces{idim} = sp_bspline (knots_hcurl{idim}, degree_hcurl{idim}, msh);
end
space = sp_vector (scalar_spaces, msh, 'curl-preserving');
sp_mul = sp_bspline (knots, degree, msh);

% Assemble the matrices
if (msh.rdim == 2)
  invmu = @(x,y) 1./c_magn_perm (x,y);
elseif (msh.rdim == 3)
  invmu = @(x,y,z) 1./c_magn_perm (x,y,z);
end

stiff_mat = op_curlu_curlv_tp (space, space, msh, invmu);
mass_mat  = op_u_v_tp (space, space, msh, c_elec_perm);
saddle_mat = op_v_gradp_tp (space, sp_mul, msh, c_elec_perm);

% Apply homogeneous Dirichlet boundary conditions
drchlt_dofs = []; drchlt_dofs_mul = [];
for iside = 1:numel (drchlt_sides)
  drchlt_dofs = union (drchlt_dofs, space.boundary(drchlt_sides(iside)).dofs);
  drchlt_dofs_mul = union (drchlt_dofs_mul, sp_mul.boundary(drchlt_sides(iside)).dofs);
end
if (isempty (drchlt_dofs_mul))
  drchlt_dofs_mul = sp_mul.ndof;
end

int_dofs = setdiff (1:space.ndof, drchlt_dofs);
int_dofs_mul = setdiff (1:sp_mul.ndof, drchlt_dofs_mul);

% Solve the eigenvalue problem
stiff_mat  = stiff_mat (int_dofs, int_dofs);
mass_mat   = mass_mat (int_dofs, int_dofs);
saddle_mat = saddle_mat (int_dofs_mul, int_dofs);

A = [stiff_mat, saddle_mat.'; ...
     saddle_mat, sparse(numel(int_dofs_mul),numel(int_dofs_mul))];
M = [mass_mat, sparse(numel(int_dofs), numel(int_dofs_mul)); ...
     sparse(numel(int_dofs_mul), numel(int_dofs)+numel(int_dofs_mul))];

eigf = zeros (space.ndof + sp_mul.ndof, numel(int_dofs) + numel(int_dofs_mul));
[eigf([int_dofs space.ndof+int_dofs_mul], :), eigv] = eig (full(A), full(M));
eigv = diag (eigv);

end

%!demo
%! ex_maxwell_eig_mixed1_square

%!demo
%! ex_maxwell_eig_mixed1_Lshaped

%!demo
%! ex_maxwell_eig_mixed1_cube

%!demo
%! ex_maxwell_eig_mixed1_thick_L
