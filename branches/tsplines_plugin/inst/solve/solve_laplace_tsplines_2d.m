% SOLVE_LAPLACE_TSPLINES_2D: Solve a 2d Laplace problem with a T-spline isoparametric discretization. 
%
% The function solves the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega = F((0,1)^2)
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
% 
% USAGE:
%
%  [tspline, msh, space, u] = solve_laplace_tsplines_2d (problem_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:  name of the file containing the T-spline Bezier extraction
%    - c_diff:    diffusion coefficient (epsilon in the equation)
%    - f:         source term
%    - g:         function for Neumann condition (optional)
%    - h:         function for Dirichlet condition
%
% OUTPUT:
%
%  tspline:  structure with the geometry and Bezier extraction (see read_bezier_extraction)
%  msh:      mesh object that defines the quadrature rule (see tspline_mesh_space)
%  space:    space object that defines the discrete space (see tspline_mesh_space)
%  u:        the computed degrees of freedom
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011, 2012, Rafael Vazquez
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

function [tspline, msh, space, u] = ...
              solve_laplace_tsplines_2d (problem_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end

% Read the T-splines plug-in file, and construct the msh and space structures
tspline = read_bezier_extraction (geo_name);
[msh, space] = tspline_mesh_space (tspline);

% Assemble the matrix and the right-hand side
x = reshape (msh.geo_map(1,:,:), msh.nqn, msh.nel);
y = reshape (msh.geo_map(2,:,:), msh.nqn, msh.nel);

stiff_mat = op_gradu_gradv (space, space, msh, c_diff(x, y));
rhs       = op_f_v (space, msh, f(x, y));

% Apply Neumann boundary conditions
for iside = nmnn_sides
  [msh_side, sp_side] = tspline_boundary (tspline, iside);

  x = reshape (msh_side.geo_map(1,:,:), msh_side.nqn, msh_side.nel);
  y = reshape (msh_side.geo_map(2,:,:), msh_side.nqn, msh_side.nel);
  gval = reshape (g (x, y, iside), msh_side.nqn, msh_side.nel);

% For T-splines the boundary fields carry the global numbers
% The field dofs is not used. This is different from other examples in GeoPDEs.
  rhs = rhs + op_f_v (sp_side, msh_side, gval);
end

% Apply Dirichlet boundary conditions
[u, drchlt_dofs] = tspline_drchlt_l2_proj (tspline, h, drchlt_sides);

int_dofs = setdiff (1:space.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u(drchlt_dofs);

% Solve the linear system
u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);

end
