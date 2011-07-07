% MP_SOLVE_STOKES_2D: Solve a Stokes flow problem on a two-dimensional multipatch domain.
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
%  [geometry, msh, space_v, vel, gnum, space_p, press, gnump] = ...
%                        mp_solve_stokes_2d (problem_data, method_data)
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
%  msh:      cell array of mesh objects (see msh_2d)
%  space_v:  cell array of space objects for the velocity (see sp_vector_2d)
%  vel:      the computed degrees of freedom for the velocity
%  gnum:     global numbering of the local degrees of freedom on each patch
%  space_p:  cell array of space objects for the pressure (see sp_bspline_2d)
%  press:    the computed degrees of freedom for the pressure
%  gnump:    global numbering of the local degrees of freedom on each patch
%             for the pressure
%
%  See also EX_STOKES_SQUARE_MP for an example
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2010, 2011 Rafael Vazquez
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
              mp_solve_stokes_2d (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

if (strcmp (lower (element_name), 'rt') || strcmp (lower (element_name), 'ndl'))
  error ('mp_solve_stokes_2d: multipatch is not ready for RT and NDL elements')
end

% Construct geometry structure, and information for interfaces and boundaries
[geometry, boundaries, interfaces] = mp_geo_load (geo_name);
npatch = numel (geometry);

ndofp = 0;
for iptc = 1:npatch
  [knotsp, knotsv1, degreev1, knotsv2, degreev2, der2] = ...
     sp_fluid_set_options_2d (element_name, geometry(iptc).nurbs.knots, ...
                              nsub, degree, regularity);

% Construct msh structure
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (knotsv1, rule);
  msh{iptc} = msh_2d (knotsv1, qn, qw, geometry(iptc), 'der2', der2);

% Construct space structure
  [spv{iptc}, spp{iptc}] = sp_bspline_fluid_2d_phys (element_name, ...
            knotsv1, degreev1, knotsv2, degreev2, knotsp, degree, msh{iptc});
end

% Create a correspondence between patches on the interfaces
[gnum,  ndof]  = mp_interface_vector_2d (interfaces, spv);
[gnump, ndofp] = mp_interface_2d (interfaces, spp);


% Compute and assemble the matrices
ncounterA = 0; ncounterM = 0; ncounterB = 0;
F = zeros (ndof, 1);

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

  [rs, cs, vs] = op_u_v_tp (spp{iptc}, spp{iptc}, msh{iptc}, @(x, y) ones (size(x))); 
  rows_M(ncounterM+(1:numel (rs))) = gnump{iptc}(rs);
  cols_M(ncounterM+(1:numel (rs))) = gnump{iptc}(cs);
  vals_M(ncounterM+(1:numel (rs))) = vs;
  ncounterM = ncounterM + numel (rs);

  F_loc = op_f_v_tp (spv{iptc}, msh{iptc}, f);
  F(gnum{iptc}) = F(gnum{iptc}) + F_loc;
end
A = sparse (rows_A, cols_A, vals_A, ndof, ndof);
B = sparse (rows_B, cols_B, vals_B, ndofp, ndof);
M = sparse (rows_M, cols_M, vals_M, ndofp, ndofp);
E = sum (M, 1) / sum (sum (M));

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
