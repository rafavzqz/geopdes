% SOLVE_STOKES_3D_BSPLINES: Solve a 3d Stokes problem with a B-spline discretization. 
%
% The function solves the stokes problem
%
%
%   -mu div grad vel + grad press = f    in Omega
%                         div vel = 0    in Omega
%                      mu dvel/dn = g    on Gamma_N
%                             vel = h    on Gamma_D
%
% USAGE:
%
%  [geometry, msh, space_v, vel, space_p, press] = ...
%                        solve_stokes_3d_bsplines (problem_data, method_data)
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
%    - viscosity:  viscosity coefficient (mu in the equation)
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:       degree of the spline functions for pressure
%    - regularity:   continuity of the spline functions for pressure
%    - n_sub:        number of subdivisions for refinement for the pressure
%    - nquad:        number of points for Gaussian quadrature rule
%    - element_name: one of {TH,SG,RT}, specify how to build the velocity space
%                    from the data for the preassure space
%                     +TH is the generalized Taylor-Hood element
%                     +SG is the SubGrid element
%                     +RT the Raviart-Thomas element
%                    for reference see ... citation
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh structure (see msh_push_forward_2d)
%  space_v:  space structure for the velocity (see sp_bspline_2d_phys)
%  vel:      the computed coeficcients of the velocity
%  space_p:  space structure for the pressure 
%  press:    the computed coeficcients of the pressure
%
% See also EX_STOKES_FLOW_3D for an example.
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


function [geometry, msh, space_v, vel, space_p, press] = ...
                          solve_stokes_3d_bsplines (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} sprintf('= problem_data.(data_names{iopt});')]);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} sprintf('= method_data.(data_names{iopt});')]);
end

% rename pressure data
degree_p     = degree;
regularity_p = regularity;
nsub_p       = nsub;

switch lower(element_name)
  case 'th'
    degree_v     = degree+1;
    regularity_v = regularity;
    nsub_v       = nsub;
  case 'sg'
    degree_v     = degree+1;
    regularity_v = regularity+1;
    nsub_v       = 2*nsub;
end


% load geometry
geometry    = geo_load (geo_name);

[knots_p, zeta_p] = kntrefine (geometry.nurbs.knots, nsub_p, degree_p, regularity_p);
[knots_v, zeta_v] = kntrefine (geometry.nurbs.knots, nsub_v, degree_v, regularity_v);

% Compute the mesh structure using the finest mesh
rule        = msh_gauss_nodes (nquad);
[qn, qw]    = msh_set_quad_nodes (zeta_v, rule);
msh         = msh_3d_tensor_product (zeta_v, qn, qw); 
msh         = msh_push_forward_3d (msh, geometry);

% Compute the space structures
space_p         = sp_bspline_3d_phys (knots_p, degree_p, msh, 'gradient', false); 
space_v         = sp_bspline_3d_phys (knots_v, degree_v, msh);
space_v         = sp_scalar_to_vector_3d (space_v, space_v, space_v, msh, 'divergence', true);

% Assemble the matrices
x           = squeeze (msh.geo_map(1,:,:));
y           = squeeze (msh.geo_map(2,:,:));
z           = squeeze (msh.geo_map(3,:,:));
%[x, y, z]   = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)), squeeze (msh.geo_map(3,:,:)));
A           = op_gradu_gradv (space_v, space_v, msh, viscosity (x, y, z)); 
B           = op_div_v_q (space_v, space_p, msh); 
M           = op_u_v (space_p, space_p, msh, ones (size (x))); 
E           = sum (M, 1) / sum (sum (M)); 
F           = op_f_v (space_v, msh, f (x, y, z)); 

vel         = zeros (space_v.ndof, 1);
press       = zeros (space_p.ndof, 1);

% Apply Dirichlet  boundary conditions
[vel_drchlt, drchlt_dofs] = sp_drchlt_l2_proj(space_v, msh, h, drchlt_sides);
vel(drchlt_dofs) = vel_drchlt;
int_dofs = setdiff (1:space_v.ndof, drchlt_dofs);
nintdofs = numel (int_dofs);
rhs_dir  = -A(int_dofs, drchlt_dofs)*vel(drchlt_dofs);

% Apply Neumann boundary conditions
rhs_nmnn = zeros(space_v.ndof,1);
for iside = nmnn_sides
  x = squeeze (msh.boundary(iside).geo_map(1,:,:));
  y = squeeze (msh.boundary(iside).geo_map(2,:,:));
  z = squeeze (msh.boundary(iside).geo_map(3,:,:));
  gval = reshape (g (x, y, z, iside), 3, msh.boundary(iside).nqn, msh.boundary(iside).nel);

  rhs_nmnn(space_v.boundary(iside).dofs) = rhs_nmnn(space_v.boundary(iside).dofs) + op_f_v(space_v.boundary(iside), msh.boundary(iside), gval);
end
rhs_nmnn = rhs_nmnn(int_dofs);

% Solve the linear system
mat = [ A(int_dofs, int_dofs), -B(:,int_dofs).',            sparse(nintdofs, 1);
       -B(:,int_dofs),          sparse(space_p.ndof, space_p.ndof), E';
        sparse(1, nintdofs),    E,                          0];
rhs = [F(int_dofs) + rhs_dir + rhs_nmnn; 
       zeros(space_p.ndof, 1); 
       0];
sol = mat \ rhs;
vel(int_dofs) = sol(1:nintdofs);
press = sol(1+nintdofs:end-1);


