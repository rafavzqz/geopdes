% SOLVE_MAXWELL_EIG: Solve the Maxwell eigenvalue problem with a B-spline discretization.
%
% Example to solve the problem
%
%    curl (1/mu(x) curl (u)) = lambda (epsilon(x) u)   in Omega = F((0,1)^n)
%      (1/mu(x) curl(u)) x n = 0                       on Gamma_N
%                      u x n = 0                       on Gamma_D
%
% USAGE:
%
%  [geometry, msh, space, eigv, eigf] = 
%                  solve_maxwell_eig (problem_data, method_data)
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
%  space:    space object that defines the discrete functions (see sp_vector)
%  eigv:     the computed eigenvalues
%  eigf:     degrees of freedom of the associated eigenfunctions
%
% See also EX_MAXWELL_EIG_SQUARE for an example
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

function [geometry, msh, space, eigv, eigf] = ...
              solve_maxwell_eig (problem_data, method_data)

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
if (msh.rdim == 2)
  invmu = @(x,y) 1./c_magn_perm (x,y);
elseif (msh.rdim == 3)
  invmu = @(x,y,z) 1./c_magn_perm (x,y,z);
end

stiff_mat = op_curlu_curlv_tp (space, space, msh, invmu);
mass_mat  = op_u_v_tp (space, space, msh, c_elec_perm);

% Apply homogeneous Dirichlet boundary conditions
drchlt_dofs = [];
for iside = 1:numel (drchlt_sides)
  drchlt_dofs = union (drchlt_dofs, space.boundary(drchlt_sides(iside)).dofs);
end
int_dofs = setdiff (1:space.ndof, drchlt_dofs);

% Solve the eigenvalue problem
eigf = zeros (space.ndof, numel(int_dofs));
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

%!demo
%! ex_maxwell_eig_cube

%!demo
%! ex_maxwell_eig_thick_L
