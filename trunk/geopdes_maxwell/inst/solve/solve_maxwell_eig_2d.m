% SOLVE_MAXWELL_EIG_2D: Solve the 2d Maxwell eigenvalue problem with a B-spline discretization.
%
% Example to solve the problem
%
%    curl (1/mu(x) curl (u)) = lambda (epsilon(x) u)   in Omega = F((0,1)^2)
%      (1/mu(x) curl(u)) x n = 0                       on Gamma_N
%                      u x n = 0                       on Gamma_D
%
% USAGE:
%
%  [geometry, msh, space, eigv, eigf] = 
%                  solve_maxwell_eig_2d (problem_data, method_data)
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
%  msh:      mesh structure (see msh_push_forward_2d)
%  space:    space structure (see sp_bspline_curl_transform_2d)
%  eigv:     the computed eigenvalues
%  eigf:     degrees of freedom of the associated eigenfunctions
%
% See also EX_MAXWELL_EIG_SQUARE for an example
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

function [geometry, msh, sp, eigv, eigf] = ...
              solve_maxwell_eig_2d (problem_data, method_data)

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

[knots, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);
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
eigv = diag (eigv);

end

%!demo
%! ex_maxwell_eig_square

%!demo
%! ex_maxwell_eig_Lshaped

%!demo
%! ex_maxwell_eig_curvedL

%!demo
%! ex_maxwell_eig_ring_1eighth


