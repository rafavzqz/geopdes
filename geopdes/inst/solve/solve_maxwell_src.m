% SOLVE_MAXWELL_SRC: Solve a Maxwell source problem with a B-spline discretization.
%
% Example to solve the problem
%
%    curl ( epsilon(x) curl (u)) + mu u = f     in Omega = F((0,1)^n)
%              (epsilon(x) curl(u)) x n = g     on Gamma_N
%                                 u x n = h x n on Gamma_D
%
% USAGE:
%
%  [geometry, msh, space, u] = solve_maxwell_src (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_mass:       coefficient for the mass matrix (mu in the equation)
%    - c_stiff:      coefficient for the stiffness matrix (epsilon in the equation)
%    - f:            source term
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
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
%  space:    space object that defines the discrete functions (see sp_vector)
%  u:        the computed degrees of freedom
%
% See also EX_MAXWELL_SRC_LSHAPED for an example.
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

function [geometry, msh, space, u] = ...
              solve_maxwell_src (problem_data, method_data)

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

% Construct space structure
scalar_spaces = cell (msh.ndim, 1);
for idim = 1:msh.ndim
  scalar_spaces{idim} = sp_bspline (knots_hcurl{idim}, degree_hcurl{idim}, msh);
end
space = sp_vector (scalar_spaces, msh, 'curl-preserving');
clear scalar_spaces

% Assemble the matrices
stiff_mat = op_curlu_curlv_tp (space, space, msh, c_stiff);
mass_mat  = op_u_v_tp (space, space, msh, c_mass);
rhs       = op_f_v_tp (space, msh, f);

% Apply Neumann boundary conditions
for iside = nmnn_sides
% Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
  gside = @(varargin) g(varargin{:},iside);
  dofs = space.boundary(iside).dofs;
  rhs(dofs) = rhs(dofs) + op_f_v_tp (space.boundary(iside), msh.boundary(iside), gside);
end

% Apply Dirichlet boundary conditions
u = zeros (space.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:space.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt ...
                              - mass_mat(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system
u(int_dofs) = (stiff_mat(int_dofs, int_dofs) + ...
                mass_mat(int_dofs, int_dofs)) \ rhs(int_dofs);

end

%!demo
%! ex_maxwell_src_square

%!demo
%! ex_maxwell_src_Lshaped

%!demo
%! ex_maxwell_src_ring

%!demo
%! ex_maxwell_src_pacman