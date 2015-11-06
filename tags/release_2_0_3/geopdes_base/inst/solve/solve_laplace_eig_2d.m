% SOLVE_LAPLACE_EIG_2D: Solve a 2d Laplace eigenproblem with a B-spline discretization (non-isoparametric approach). 
%
% The function solves the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = lambda (mu(x) u)    in Omega = F((0,1)^2)
%                epsilon(x) du/dn = 0    on Gamma_N
%                               u = 0    on Gamma_D
%
% USAGE:
%
%  [geometry, msh, space, lambda, u] = solve_laplace_eig_2d (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_diff:       diffusion coefficient (epsilon in the equation)
%    - c_mass:       mass coefficient (mu in the equation)
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
%  space:    space object that defines the discrete space (see sp_bspline_2d)
%  lambda:   the computed eigenvalues
%  u:        degrees of freedom of the computed eigenvectors
%
% See also EX_LAPLACE_EIG_SQUARE for an example.
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011, Rafael Vazquez
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

function [geometry, msh, space, lambda, u] = ...
              solve_laplace_eig_2d (problem_data, method_data)

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
geometry  = geo_load (geo_name);

[knots, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);
  
% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_2d (zeta, qn, qw, geometry);
  
% Construct space structure
space    = sp_bspline_2d (knots, degree, msh);
  
% Assemble the matrices
stiff_mat = op_gradu_gradv_tp (space, space, msh, c_diff);
mass_mat  = op_u_v_tp (space, space, msh, c_mass);

% Apply homogeneous Dirichlet boundary conditions
drchlt_dofs = [];
for iside = 1:numel (drchlt_sides)
  drchlt_dofs = union (drchlt_dofs, space.boundary(drchlt_sides(iside)).dofs);
end
int_dofs = setdiff (1:space.ndof, drchlt_dofs);

% Solve the eigenvalue problem
u = zeros (space.ndof, numel(int_dofs));
[u(int_dofs, :), lambda] = ...
    eig (full (stiff_mat(int_dofs, int_dofs)), full (mass_mat(int_dofs, int_dofs)));
lambda = diag (lambda);

end

%!demo
%! ex_laplace_eig_square
