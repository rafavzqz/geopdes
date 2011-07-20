% MP_SOLVE_STOKES_3D: Solve a Stokes flow problem on a three-dimensional multipatch domain.
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
%                        mp_solve_stokes_3d (problem_data, method_data)
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
%  msh:      array of mesh structures (see msh_push_forward_3d)
%  space_v:  array of space structures for the velocity (see sp_bspline_3d_phys)
%  vel:      the computed coeficcients of the velocity
%  gnum:     global numbering of the local degrees of freedom on each patch
%  space_p:  array of space structures for the pressure 
%  press:    the computed coeficcients of the pressure
%  gnump:    global numbering of the local degrees of freedom on each patch
%             for the pressure
%
%  See also EX_STOKES_DRIVEN_CAVITY_3D_MP for an example
%
% In order to make it work, a script file with the data must also be executed.
% The available data files for this problem are:
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
              mp_solve_stokes_3d (problem_data, method_data)

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
  error ('mp_solve_stokes_3d: multipatch is not ready for RT and NDL elements')
end

% Construct geometry structure, and information for interfaces and boundaries
[geometry, boundaries, interfaces] = mp_geo_load (geo_name);
npatch = numel (geometry);

ndofp = 0;
for iptc = 1:npatch
  [knotsp, knotsv1, degreev1, knotsv2, degreev2, knotsv3, degreev3, der2] = ...
     sp_fluid_set_options_3d (element_name, geometry(iptc).nurbs.knots, ...
                              nsub, degree, regularity);

% Construct msh structure
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (knotsv1, rule);
  msh{iptc} = msh_3d_tensor_product (knotsv1, qn, qw);
  msh{iptc} = msh_push_forward_3d (msh{iptc}, geometry(iptc));

% Construct space structure
  [spv{iptc}, spp{iptc}] = sp_bspline_fluid_3d_phys (element_name, ...
            knotsv1, degreev1, knotsv2, degreev2, knotsv3, degreev3, ...
            knotsp, degree, msh{iptc});
end

% Create a correspondence between patches on the interfaces
[gnum,  ndof]  = mp_interface_vector_3d (interfaces, spv);
[gnump, ndofp] = mp_interface_3d (interfaces, spp); 

% Compute and assemble the matrices
nent = sum (cellfun (@(x, y, z) x.nel * y.nsh_max * z.nsh_max, msh, spv, spv));
rows_A = zeros (nent, 1); cols_A = zeros (nent, 1); vals_A = zeros (nent, 1);
ncounterA = 0;

nent = sum (cellfun (@(x, y, z) x.nel * y.nsh_max * z.nsh_max, msh, spp, spp));
rows_M = zeros (nent, 1); cols_M = zeros (nent, 1); vals_M = zeros (nent, 1);
ncounterM = 0;

nent = sum (cellfun (@(x, y, z) x.nel * y.nsh_max * z.nsh_max, msh, spv, spp));
rows_B = zeros (nent, 1); cols_B = zeros (nent, 1); vals_B = zeros (nent, 1);
ncounterB = 0;

F = zeros (ndof, 1);

for iptc = 1:npatch
  x = squeeze (msh{iptc}.geo_map(1,:,:));
  y = squeeze (msh{iptc}.geo_map(2,:,:));
  z = squeeze (msh{iptc}.geo_map(3,:,:));
  visc = reshape (viscosity (x, y, z), msh{iptc}.nqn, msh{iptc}.nel);
  fval = reshape (f (x, y, z), 3, msh{iptc}.nqn, msh{iptc}.nel);

  [rs, cs, vs] = op_gradu_gradv (spv{iptc}, spv{iptc}, msh{iptc}, visc); 
  rows_A(ncounterA+(1:numel (rs))) = gnum{iptc}(rs);
  cols_A(ncounterA+(1:numel (rs))) = gnum{iptc}(cs);
  vals_A(ncounterA+(1:numel (rs))) = vs;
  ncounterA = ncounterA + numel (rs);

  [rs, cs, vs] = op_div_v_q (spv{iptc}, spp{iptc}, msh{iptc}); 
  rows_B(ncounterB+(1:numel (rs))) = gnump{iptc}(rs);
  cols_B(ncounterB+(1:numel (rs))) = gnum{iptc}(cs);
  vals_B(ncounterB+(1:numel (rs))) = vs;
  ncounterB = ncounterB + numel (rs);

  [rs, cs, vs] = op_u_v (spp{iptc}, spp{iptc}, msh{iptc}, ones (msh{iptc}.nqn, msh{iptc}.nel)); 
  rows_M(ncounterM+(1:numel (rs))) = gnump{iptc}(rs);
  cols_M(ncounterM+(1:numel (rs))) = gnump{iptc}(cs);
  vals_M(ncounterM+(1:numel (rs))) = vs;
  ncounterM = ncounterM + numel (rs);

  F_loc = op_f_v (spv{iptc}, msh{iptc}, fval);
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