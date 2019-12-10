% SOLVE_KIRCHHOFF_LOVE_SHELL: Solve the Kirchhoff-Love shell model in a NURBS domain.
%
% USAGE:
%
%  [geometry, msh, space, u] = solve_kirchhoff_love_shell (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - drchlt_components: cell-array, the components that are set to zero for each drchlt_side
%    - E_coeff:      function handle for Young's modulus
%    - nu_coeff:     function handle for Poisson's ratio
%    - thickness:    scalar value, thickness of the shell
%    - f:            source term, distributed load
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree (>=2) of the spline functions.
%    - regularity: continuity (>=1) of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_cartesian)
%  space:    space object that defines the discrete basis functions (see sp_vector)
%  u:        the computed degrees of freedom
%
% NOTE: only homogeneous Dirichlet conditions implemented so far
%
% See also EX_KL_SHELL_SCORDELIS_LO_ROOF for an example.
%
% Copyright (C) 2017-2019 Pablo Antolin, Luca Coradello, Rafael Vazquez
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

function [geometry, msh, space, u] = solve_kirchhoff_love_shell (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

if (any(degree <= 1) || any(regularity == 0))
  error ('The degree must be at least two, and the regularity at least C^1')
end

geometry = geo_load (geo_name);
degelev  = max (degree - (geometry.nurbs.order-1), 0);
nurbs    = nrbdegelev (geometry.nurbs, degelev);
[rknots, ~, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);

nurbs    = nrbkntins (nurbs, nknots);
geometry = geo_load (nurbs);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (geometry.nurbs.knots, rule);
msh      = msh_cartesian (geometry.nurbs.knots, qn, qw, geometry);

% Construct space structure
sp_scalar = sp_nurbs (geometry.nurbs, msh);
scalar_spaces = repmat ({sp_scalar}, 1, msh.rdim);
space = sp_vector (scalar_spaces, msh);

% Assemble the stiffness matrix and right-hand side
K = op_KL_shells_tp (space, space, msh, E_coeff, nu_coeff, thickness);
rhs = op_f_v_tp (space, msh, f);

u = zeros (space.ndof, 1);

% Apply boundary conditions
drchlt_dofs = [];
for iside = 1:numel(drchlt_sides)
  side = drchlt_sides(iside);
  if (~exist('drchlt_components','var'))
    components = 1:3;
  else
    components = drchlt_components{iside};
  end
  for icomp = components
    drchlt_dofs = union (drchlt_dofs, space.boundary(side).dofs(space.boundary(side).comp_dofs{icomp}));
  end
end

int_dofs = setdiff (1:space.ndof, drchlt_dofs);
%rhs(int_dofs) = rhs(int_dofs) - K(int_dofs, drchlt_dofs)*u(drchlt_dofs);

% Solve the linear system
u(int_dofs) = K(int_dofs, int_dofs) \ rhs(int_dofs);

end
