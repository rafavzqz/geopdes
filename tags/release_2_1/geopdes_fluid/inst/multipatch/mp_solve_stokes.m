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
%  [geometry, msh, space_v, vel, gnum, space_p, press, gnump] = ...
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
%  msh:      cell array of mesh objects (see msh_cartesian)
%  space_v:  cell array of space objects for the velocity (see sp_vector or sp_vector_div_transform)
%  vel:      the computed degrees of freedom for the velocity
%  gnum:     global numbering of the local degrees of freedom on each patch
%  space_p:  cell array of space objects for the pressure (see sp_bspline)
%  press:    the computed degrees of freedom for the pressure
%  gnump:    global numbering of the local degrees of freedom on each patch
%             for the pressure
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

function [geometry, msh, spv, vel, gnum, spp, press, gnump] = ...
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
[geometry, boundaries, interfaces] = mp_geo_load (geo_name);
npatch = numel (geometry);

for iptc = 1:npatch
% Construct msh structure using the finest mesh
  msh_breaks = msh_set_breaks (element_name, ...
                                       geometry(iptc).nurbs.knots, nsub);
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (msh_breaks, rule);
  msh{iptc} = msh_cartesian (msh_breaks, qn, qw, geometry(iptc));

% Construct space structure
  [spv{iptc}, spp{iptc}] = sp_bspline_fluid (element_name, ...
               geometry(iptc).nurbs.knots, nsub, degree, regularity, msh{iptc});
end

% Create a correspondence between patches on the interfaces
[gnum,  ndof]  = mp_interface_vector (interfaces, spv);
[gnump, ndofp] = mp_interface (interfaces, spp);

% Compute and assemble the matrices
ncounterA = 0; ncounterM = 0; ncounterB = 0;
E = zeros (1, ndofp); F = zeros (ndof, 1);

if (msh{1}.rdim == 2)
  fun_one = @(x, y) ones (size(x));
elseif (msh{1}.rdim == 3)
  fun_one = @(x, y, z) ones (size(x));
end

for iptc = 1:npatch
  [rs, cs, vs] = op_gradu_gradv_tp (spv{iptc}, spv{iptc}, msh{iptc}, viscosity); 
  rows_A(ncounterA+(1:numel (rs))) = gnum{iptc}(rs);
  cols_A(ncounterA+(1:numel (rs))) = gnum{iptc}(cs);
  vals_A(ncounterA+(1:numel (rs))) = vs;
  ncounterA = ncounterA + numel (rs);

  [rs, cs, vs] = op_div_v_q_tp (spv{iptc}, spp{iptc}, msh{iptc}); 
  rows_B(ncounterB+(1:numel (rs))) = gnump{iptc}(rs);
  cols_B(ncounterB+(1:numel (rs))) = gnum{iptc}(cs);
  vals_B(ncounterB+(1:numel (rs))) = vs;
  ncounterB = ncounterB + numel (rs);

  E_loc = op_f_v_tp (spp{iptc}, msh{iptc}, fun_one).';
  E(gnump{iptc}) = E(gnump{iptc}) + E_loc;

  F_loc = op_f_v_tp (spv{iptc}, msh{iptc}, f);
  F(gnum{iptc}) = F(gnum{iptc}) + F_loc;
end

clear rs cs vs
A = sparse (rows_A, cols_A, vals_A, ndof, ndof);
clear rows_A cols_A vals_A
B = sparse (rows_B, cols_B, vals_B, ndofp, ndof);
clear rows_B cols_B vals_B

vel   = zeros (ndof, 1);
press = zeros (ndofp, 1);

% Apply Dirichlet boundary conditions
[vel_drchlt, drchlt_dofs] = mp_sp_drchlt_l2_proj (spv, msh, h, gnum, boundaries, drchlt_sides);
vel(drchlt_dofs) = vel_drchlt;

int_dofs = setdiff (1:ndof, drchlt_dofs);
nintdofs = numel (int_dofs);

% Solve the linear system
mat = [A(int_dofs, int_dofs), -B(:,int_dofs).', sparse(nintdofs, 1);
       -B(:,int_dofs), sparse(ndofp, ndofp), E';
       sparse(1, nintdofs), E, 0];
rhs = [F(int_dofs)-A(int_dofs, drchlt_dofs)*vel(drchlt_dofs); 
       B(:, drchlt_dofs)*vel(drchlt_dofs); 
       0];

sol = mat \ rhs;

vel(int_dofs) = sol(1:nintdofs);
press = sol(1+nintdofs:end-1);

end
