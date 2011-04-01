% SOLVE_STOKES_2D_BSPLINES: Solve a 2d Stokes problem with a B-spline discretization. 
%
% The function solves the Stokes problem
%
%   -div(mu(x) grad(vel)) + grad(press) = f    in Omega
%                              div(vel) = 0    in Omega
%                         mu(x) dvel/dn = g    on Gamma_N
%                                   vel = h    on Gamma_D
%
% USAGE:
%
%  [geometry, msh, space_v, vel, space_p, press] = ...
%                        solve_stokes_2d_bsplines (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - f:            force term
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
%    - viscosity:    viscosity coefficient (mu in the equation)
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:       degree of the spline functions for pressure
%    - regularity:   continuity of the spline functions for pressure
%    - nsub:         number of subdivisions for refinement for the pressure
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
%  msh:      mesh structure (see msh_push_forward_2d)
%  space_v:  space structure for the velocity (see sp_bspline_2d_phys)
%  vel:      the computed coeficcients of the velocity
%  space_p:  space structure for the pressure 
%  press:    the computed coeficcients of the pressure
%
% See also EX_STOKES_SQUARE_* EX_STOKES_ANNULUS_* for an examples.
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011, Andrea Bressan, Rafael Vazquez
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
                          solve_stokes_2d_bsplines (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} sprintf('= problem_data.(data_names{iopt});')]);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} sprintf('= method_data.(data_names{iopt});')]);
end

% load geometry
geometry    = geo_load (geo_name);
[knotsp, knotsv1, degreev1, knotsv2, degreev2, der2] = ...
   sp_fluid_set_options_2d (element_name, geometry.nurbs.knots, nsub, degree, regularity);

% Compute the mesh structure using the finest mesh
rule        = msh_gauss_nodes (nquad);
[qn, qw]    = msh_set_quad_nodes (knotsv1, rule);
msh         = msh_2d_tensor_product (knotsv1, qn, qw);
msh         = msh_push_forward_2d (msh, geometry, 'der2', der2);

% Compute the space structures
[space_v, space_p, PI] = sp_bspline_fluid_2d_phys (element_name, ...
          knotsv1, degreev1, knotsv2, degreev2, knotsp, degree, msh);

% Assemble the matrices
x           = squeeze (msh.geo_map(1,:,:));
y           = squeeze (msh.geo_map(2,:,:));

A           = op_gradu_gradv (space_v, space_v, msh, viscosity (x, y)); 
B           = PI' * op_div_v_q (space_v, space_p, msh);
M           = op_u_v (space_p, space_p, msh, ones (size (x))); 
E           = sum (M, 1) * PI / sum (sum (M)); 
F           = op_f_v (space_v, msh, f (x, y)); 

vel         = zeros (space_v.ndof, 1);
press       = zeros (space_p.ndof, 1);

% Apply Dirichlet  boundary conditions
[vel_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space_v, msh, h, drchlt_sides);
vel(drchlt_dofs) = vel_drchlt;
int_dofs = setdiff (1:space_v.ndof, drchlt_dofs);
nintdofs = numel (int_dofs);
rhs_dir  = -A(int_dofs, drchlt_dofs)*vel(drchlt_dofs);

% Apply Neumann boundary conditions
rhs_nmnn = zeros(space_v.ndof,1);
for iside = nmnn_sides
  x = squeeze (msh.boundary(iside).geo_map(1,:,:));
  y = squeeze (msh.boundary(iside).geo_map(2,:,:));
  gval = reshape (g (x, y, iside), 2, msh.boundary(iside).nqn, msh.boundary(iside).nel);

  rhs_nmnn(space_v.boundary(iside).dofs) = rhs_nmnn(space_v.boundary(iside).dofs) + op_f_v(space_v.boundary(iside), msh.boundary(iside), gval);
end
rhs_nmnn = rhs_nmnn(int_dofs);

% Solve the linear system
mat = [ A(int_dofs, int_dofs), -B(:,int_dofs).',            sparse(nintdofs, 1);
       -B(:,int_dofs),          sparse(size (B,1), size(B,1)), E';
        sparse(1, nintdofs),    E,                          0];
rhs = [F(int_dofs) + rhs_dir + rhs_nmnn; 
       B(:, drchlt_dofs)*vel(drchlt_dofs);
       0];
sol = mat \ rhs;
vel(int_dofs) = sol(1:nintdofs);
press = PI * sol(1+nintdofs:end-1);

end
 