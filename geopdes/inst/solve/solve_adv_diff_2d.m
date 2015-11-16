% SOLVE_ADV_DIFF_2D: Solve a 2d advection-diffusion problem with a B-spline discretization. 
% The function offers the possibility to use the SUPG method (stab=true/false).
% Dirichlet boundary conditions 
%
% The function solve the advection-diffusion problem
%
%    - div (mu(x) grad (u)) + div ( vel * u ) = f  in Omega
%                                 mu(x) du/dn = g  on Gamma_N
%                                           u = h  on Gamma_D
%	
% USAGE:
%
%  [geometry, msh, space, u] = solve_adv_diff_2d (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%  - geo_name:     name of the file containing the geometry
%  - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%  - drchlt_sides: sides with Dirichlet boundary condition
%  - c_diff:       diffusion coefficient (mu in the equation)
%  - f:            source term
%  - g:            function for Neumann condition (if nmnn_sides is not empty)
%  - h:            function for Dirichlet boundary condition
%
%  method_data : a structure with discretization data. Its fields are:
%  - degree:     degree of the spline functions.
%  - regularity: continuity of the spline functions.
%  - nsub:       number of subelements with respect to the geometry mesh 
%                 (nsub=1 leaves the mesh unchanged)
%  - nquad:      number of points for Gaussian quadrature rule
%  - stab:       logical variable to decide whether to use SUPG stabilization
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_cartesian)
%  space:    space object that defines the discrete space (see sp_scalar)
%  u:        the computed degrees of freedom
%
% See also EX_ADVECTION_DIFFUSION_SQUARE for an example.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% C. De Falco, A. Reali, R. Vazquez
% GeoPDEs: a research tool for IsoGeometric Analysis of PDEs
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2014 Rafael Vazquez
% Copyright (C) 2013 Anna Tagliabue
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
              solve_adv_diff_2d (problem_data, method_data)

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
msh      = msh_cartesian (zeta, qn, qw, geometry, 'der2', true);

% Construct space structure
space    = sp_bspline (knots, degree, msh);

% Assemble the matrix and the right-hand side
mat = op_gradu_gradv_tp (space, space, msh, c_diff);
mat = mat + op_vel_dot_gradu_v_tp (space, space, msh, vel);
rhs = op_f_v_tp (space, msh, f);

% Add stabilization terms to the matrix and the right-hand-side
if (stab)
  mat = mat + op_mat_stab_SUPG_tp (space, space, msh, c_diff, grad_diff, vel);
  rhs = rhs + op_rhs_stab_SUPG_tp (space, msh, c_diff, vel, f);
end

% Apply Neumann boundary conditions
if (exist ('nmnn_sides', 'var'))
  for iside = nmnn_sides
% Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
    gside = @(varargin) g(varargin{:},iside);
    dofs = space.boundary(iside).dofs;
    rhs(dofs) = rhs(dofs) + op_f_v_tp (space.boundary(iside), msh.boundary(iside), gside);
  end
end

% Apply Dirichlet boundary conditions
u = zeros (space.ndof, 1);

[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:space.ndof, drchlt_dofs);

% % Imposing the boundary conditions without the L2 projection
% u(space.boundary(4).dofs) = 1;
% u(space.boundary(1).dofs) = 0;
% u(space.boundary(3).dofs) = 0;
% dofs = space.boundary(2).dofs;
% [~, ind] = max (u(dofs));
% u(dofs(1:ind-1)) = 0;
% u(dofs(ind:end)) = 1;

rhs(int_dofs) = rhs(int_dofs) - mat(int_dofs, drchlt_dofs)*u(drchlt_dofs);


% Solve the linear system
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

end

%!demo
%! ex_advection_diffusion_square
