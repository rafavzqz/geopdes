% SOLVE_BILAPLACE_2D_ISO: Solve a 2d bilaplace problem with a variational formulation based on the Laplacian.
%
% The function solves the bilaplacian problem
%
%   laplace( epsilon(x) laplace(u)) = f    in Omega = F((0,1)^2)
%                             du/dn = 0    on Gamma_D
%                                 u = 0    on Gamma_D
%
% with a variational formulation based on the Laplacian, that is,
%
%       (epsilon laplace(u), laplace(v)) = (f,v)
%
% USAGE:
%
%  [geometry, msh, space, u] = solve_bilaplace_2d_iso (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - drchlt_sides: sides with Dirichlet boundary condition (always homogeneous)
%    - c_diff:       physical parameter (epsilon in the equation)
%    - f:            source term
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
%  space:    space object that defines the discrete space (see sp_scalar)
%  u:        the computed degrees of freedom
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011, 2013, 2015 Rafael Vazquez
% Copyright (C) 2013, Marco Pingaro
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
              solve_bilaplace_2d_iso (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

if (any (regularity < 1))
  error ('The regularity should be at least 1')
end

% Construct geometry structure
geometry = geo_load (geo_name);
degelev  = max (degree - (geometry.nurbs.order-1), 0);
nurbs    = nrbdegelev (geometry.nurbs, degelev);
[rknots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);

nurbs = nrbkntins (nurbs, nknots);
geometry = geo_load (nurbs);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry,'der2', true);
  
% Construct space structure
space = sp_nurbs (geometry.nurbs, msh);

% Assemble the matrices
stiff_mat = op_laplaceu_laplacev_tp (space, space, msh, c_diff);
rhs       = op_f_v_tp (space, msh, f);

% Apply boundary conditions
% Only homogeneous conditions are implemented
u = zeros (space.ndof, 1);

drchlt_dofs_u = []; drchlt_dofs_r = [];
for iside = drchlt_sides
  drchlt_dofs_u = union (drchlt_dofs_u, space.boundary(iside).dofs);
  drchlt_dofs_r = union (drchlt_dofs_r, space.boundary(iside).adjacent_dofs);
end
drchlt_dofs = union (drchlt_dofs_u, drchlt_dofs_r);

int_dofs = setdiff (1:space.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u(drchlt_dofs);

% Solve the linear system
u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);

end
