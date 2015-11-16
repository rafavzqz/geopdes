% SOLVE_STOKES: Solve a Stokes problem with a B-spline discretization. 
%
% The function solves the Stokes problem
%
%   -div(mu(x) grad(vel)) + grad(press) = f    in Omega
%                              div(vel) = 0    in Omega
%             mu(x) dvel/dn - press * n = g    on Gamma_N
%                                   vel = h    on Gamma_D
%
% USAGE:
%
%  [geometry, msh, space_v, vel, space_p, press] = ...
%                        solve_stokes (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with natural boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - f:            force term
%    - g:            function for natural condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
%    - viscosity:    viscosity coefficient (mu in the equation)
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:       degree of the spline functions for pressure
%    - regularity:   continuity of the spline functions for pressure
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   for the pressure space (nsub=1 leaves the mesh unchanged)
%    - nquad:        number of points for Gaussian quadrature rule
%    - element_name: one of {TH,SG,RT,NDL}, specify how to build the velocity
%                    space from the data for the pressure space
%                     +TH  is the generalized Taylor-Hood element
%                     +SG  is the SubGrid element
%                     +RT  the generalized Raviart-Thomas element
%                     +NDL the generalized Nedelec element of the second class
%                    for more details see the references below
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_cartesian)
%  space_v:  space object for the velocity (see sp_vector)
%  vel:      the computed degrees of freedom for the velocity
%  space_p:  space object for the pressure (see sp_scalar)
%  press:    the computed degrees of freedom for the pressure
%
% See also EX_STOKES_SQUARE_* EX_STOKES_ANNULUS_* for some examples.
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011, Andrea Bressan, Rafael Vazquez
% Copyright (C) 2014, 2015  Rafael Vazquez
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

function [geometry, msh, space_v, vel, space_p, press] = ...
                          solve_stokes (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% load geometry
geometry = geo_load (geo_name);

% Compute the mesh structure using the finest mesh
switch (upper(element_name))
  case {'RT', 'TH', 'NDL'}
    [~, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);
  case {'SG'}
    [~, zeta] = kntrefine (geometry.nurbs.knots, 2*nsub-1, degree, regularity);
end
rule       = msh_gauss_nodes (nquad);
[qn, qw]   = msh_set_quad_nodes (zeta, rule);
msh        = msh_cartesian (zeta, qn, qw, geometry);

% Compute the space structures
[space_v, space_p] = sp_bspline_fluid (element_name, ...
                geometry.nurbs.knots, nsub, degree, regularity, msh);

% Assemble the matrices
if (msh.rdim == 2)
  fun_one = @(x, y) ones (size(x));
elseif (msh.rdim == 3)
  fun_one = @(x, y, z) ones (size(x));
end
A = op_gradu_gradv_tp (space_v, space_v, msh, viscosity); 
B = op_div_v_q_tp (space_v, space_p, msh);
E = op_f_v_tp (space_p, msh, fun_one).';
F = op_f_v_tp (space_v, msh, f);

vel   = zeros (space_v.ndof, 1);
press = zeros (space_p.ndof, 1);

% Apply natural boundary conditions
rhs_nmnn = zeros(space_v.ndof,1);
for iside = nmnn_sides
  msh_side = msh_eval_boundary_side (msh, iside);
  if (strcmpi (element_name, 'RT') || strcmpi (element_name, 'NDL'))
    msh_side_from_interior = msh_boundary_side_from_interior (msh, iside);

    sp_side = space_v.constructor (msh_side_from_interior);
    sp_side = struct (sp_precompute (sp_side, msh_side_from_interior, 'value', true));
    sp_side.dofs = 1:sp_side.ndof;
  else
    sp_side  = sp_eval_boundary_side (space_v, msh_side);
  end

  x = cell (msh_side.rdim, 1);
  for idim = 1:msh_side.rdim
    x{idim} = squeeze (msh_side.geo_map(idim,:,:));
  end
  gval = reshape (g (x{:}, iside), msh.rdim, msh_side.nqn, msh_side.nel);

  rhs_nmnn(sp_side.dofs) = rhs_nmnn(sp_side.dofs) + op_f_v (sp_side, msh_side, gval);
end

% Apply Dirichlet  boundary conditions. For RT and NDL elements the normal
%  component is imposed strongly, and the tangential one is imposed weakly.
if (strcmpi (element_name, 'RT') || strcmpi (element_name, 'NDL'))
  [N_mat, N_rhs] = ...
    sp_weak_drchlt_bc (space_v, msh, drchlt_sides, h, viscosity, Cpen);
  [vel_drchlt, drchlt_dofs] = sp_drchlt_l2_proj_udotn (space_v, msh, drchlt_sides, h);
else
  [vel_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space_v, msh, h, drchlt_sides);
  N_mat = sparse (space_v.ndof, space_v.ndof, 1);
  N_rhs = zeros (space_v.ndof, 1);
end

vel(drchlt_dofs) = vel_drchlt;
int_dofs = setdiff (1:space_v.ndof, drchlt_dofs);
nintdofs = numel (int_dofs);
rhs_dir  = -A(int_dofs, drchlt_dofs)*vel(drchlt_dofs) + N_mat(int_dofs,drchlt_dofs)*vel(drchlt_dofs);

% Solve the linear system
if (isempty (nmnn_sides))
% If all the sides are Dirichlet, the pressure is zero averaged.
  mat = [ A(int_dofs, int_dofs) - N_mat(int_dofs, int_dofs), -B(:,int_dofs).',              sparse(nintdofs, 1);
         -B(:,int_dofs),                                     sparse(size (B,1), size(B,1)), E.';
          sparse(1, nintdofs),                               E,                             0];
  rhs = [F(int_dofs) + N_rhs(int_dofs) + rhs_dir; 
         B(:, drchlt_dofs)*vel(drchlt_dofs);
         0];

  sol = mat \ rhs;
  vel(int_dofs) = sol(1:nintdofs);
  press = sol(1+nintdofs:end-1);
else
% With natural boundary condition, the constraint on the pressure is not needed.
  mat = [ A(int_dofs, int_dofs) - N_mat(int_dofs, int_dofs), -B(:,int_dofs).';
         -B(:,int_dofs),                                     sparse(size (B,1), size (B,1))];
  rhs = [F(int_dofs) + N_rhs(int_dofs) + rhs_dir + rhs_nmnn(int_dofs);
         B(:, drchlt_dofs)*vel(drchlt_dofs)];

  sol = mat \ rhs;
  vel(int_dofs) = sol(1:nintdofs);
  press = sol(1+nintdofs:end);
end

end
