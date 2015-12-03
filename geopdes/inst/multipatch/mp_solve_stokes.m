% MP_SOLVE_STOKES: Solve a Stokes flow problem on a multipatch domain.
%
% The function solves the Stokes problem
%
%   -div(mu(x) grad(vel)) + grad(press) = f    in Omega
%                              div(vel) = 0    in Omega
%                                   vel = h    on Gamma_D
%
% USAGE:
%
%  [geometry, msh, space_v, vel, space_p, press] = ...
%                        mp_solve_stokes (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - f:            force term
%    - h:            function for Dirichlet boundary condition
%    - viscosity:    viscosity coefficient (mu in the equation)
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:       degree of the spline functions for pressure
%    - regularity:   continuity of the spline functions for pressure
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   for the pressure space (nsub=1 leaves the mesh unchanged)
%    - nquad:        number of points for Gaussian quadrature rule
%    - element_name: one of {TH,SG}, specify how to build the velocity
%                    space from the data for the pressure space
%                     +TH  is the generalized Taylor-Hood element
%                     +SG  is the SubGrid element
%
% OUTPUT:
%
%  geometry: array of geometry structures (see geo_load)
%  msh:      multipatch mesh, consisting of several Cartesian meshes (see msh_multipatch)
%  space_v:  multipatch space, formed by several tensor product spaces plus the connectivity (see sp_multipatch)
%  vel:      the computed degrees of freedom for the velocity
%  space_p:  multipatch space for the pressure (see sp_multipatch)
%  press:    the computed degrees of freedom for the pressure
%
%  See also EX_STOKES_DRIVEN_CAVITY_3D_MP for an example
%
% Copyright (C) 2009, 2010 Carlo de Falco
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

function [geometry, msh, space_v, vel, space_p, press] = ...
              mp_solve_stokes (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

if (strcmpi ('element_name', 'RT') || strcmpi ('element_name', 'NDL'))
  error ('For RT and NDL spaces, use mp_solve_stokes_div_conforming');
end

% Construct geometry structure, and information for interfaces and boundaries
[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (geo_name);
npatch = numel (geometry);

msh = cell (1, npatch);
spv = cell (1, npatch);
spp = cell (1, npatch);
for iptc = 1:npatch
% Construct msh structure using the finest mesh
  msh_breaks = msh_set_breaks (element_name, geometry(iptc).nurbs.knots, nsub);
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (msh_breaks, rule);
  msh{iptc} = msh_cartesian (msh_breaks, qn, qw, geometry(iptc));

% Construct space structure
  [spv{iptc}, spp{iptc}] = sp_bspline_fluid (element_name, ...
               geometry(iptc).nurbs.knots, nsub, degree, regularity, msh{iptc});
end

msh = msh_multipatch (msh, boundaries);
space_v = sp_multipatch (spv, msh, interfaces, boundary_interfaces);
space_p = sp_multipatch (spp, msh, interfaces, boundary_interfaces);
clear spv spp

% Compute and assemble the matrices
if (msh.rdim == 2)
  fun_one = @(x, y) ones (size(x));
elseif (msh.rdim == 3)
  fun_one = @(x, y, z) ones (size(x));
end

A = op_gradu_gradv_mp (space_v, space_v, msh, viscosity);
B = op_div_v_q_mp (space_v, space_p, msh);
E = (op_f_v_mp (space_p, msh, fun_one)).';
F = op_f_v_mp (space_v, msh, f);

vel   = zeros (space_v.ndof, 1);
press = zeros (space_p.ndof, 1);

% Apply Dirichlet boundary conditions
[vel_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space_v, msh, h, drchlt_sides);
vel(drchlt_dofs) = vel_drchlt;

int_dofs = setdiff (1:space_v.ndof, drchlt_dofs);
nintdofs = numel (int_dofs);

% Solve the linear system
mat = [A(int_dofs, int_dofs), -B(:,int_dofs).', sparse(nintdofs, 1);
       -B(:,int_dofs), sparse(space_p.ndof, space_p.ndof), E.';
       sparse(1, nintdofs), E, 0];
rhs = [F(int_dofs)-A(int_dofs, drchlt_dofs)*vel(drchlt_dofs); 
       B(:, drchlt_dofs)*vel(drchlt_dofs); 
       0];

sol = mat \ rhs;

vel(int_dofs) = sol(1:nintdofs);
press = sol(1+nintdofs:end-1);

end
