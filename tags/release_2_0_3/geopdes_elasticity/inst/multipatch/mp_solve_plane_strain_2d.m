% MP_SOLVE_PLANE_STRAIN_2D: Solve a plane-strain problem on a two-dimensional multipatch domain.
%
% Example to solve the plane strain problem
%
%      - div (sigma(u)) = f    in Omega
%      sigma(u) \cdot n = g    on Gamma_N
%                     u = h    on Gamma_D
%
% with   sigma(u) = mu*(grad(u) + grad(u)^t) + lambda*div(u)*I,
% and the domain \Omega is formed by several patches of the form F((0,1)^2).
%
%   u:          displacement vector
%   sigma:      Cauchy stress tensor
%   lambda, mu: Lame' parameters
%   I:          identity tensor
%
% USAGE:
%
%  [geometry, msh, space, u, gnum] = 
%             mp_solve_plane_strain_2d (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - lambda_lame:  first Lame' parameter
%    - mu_lame:      second Lame' parameter
%    - f:            source term
%    - h:            function for Dirichlet boundary condition
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
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
%  geometry: array of geometry structures (see geo_load)
%  msh:      cell array of mesh objects (see msh_2d)
%  space:    cell array of space object (see sp_vector_2d)
%  u:        the computed degrees of freedom
%  gnum:     global numbering of the degrees of freedom
%
% Copyright (C) 2010, 2011 Carlo de Falco, Rafael Vazquez
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

function [geometry, msh, sp, u, gnum] = ...
              mp_solve_plane_strain_2d (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Construct geometry structure, and information for interfaces and boundaries
[geometry, boundaries, interfaces] = mp_geo_load (geo_name);
npatch = numel (geometry);

for iptc = 1:npatch
  degelev  = max (degree - (geometry(iptc).nurbs.order-1), 0);
  nurbs    = nrbdegelev (geometry(iptc).nurbs, degelev);
  [rknots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);

  nurbs    = nrbkntins (nurbs, nknots);
  geometry(iptc) = orderfields (geo_load (nurbs), geometry);

% Construct msh structure
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (geometry(iptc).nurbs.knots, rule);
  msh{iptc} = msh_2d (geometry(iptc).nurbs.knots, qn, qw, geometry(iptc));

% Construct space structure
  sp_scalar = sp_nurbs_2d (geometry(iptc).nurbs, msh{iptc});
  sp{iptc} = sp_vector_2d (sp_scalar, sp_scalar, msh{iptc});
end
clear sp_scalar

% Create a correspondence between patches on the interfaces
[gnum, ndof] = mp_interface_vector_2d (interfaces, sp);

% Compute and assemble the matrices
rhs = zeros (ndof, 1);

ncounter = 0;
for iptc = 1:npatch
  [rs, cs, vs] = op_su_ev_tp (sp{iptc}, sp{iptc}, msh{iptc}, lambda_lame, mu_lame);

  rows(ncounter+(1:numel (rs))) = gnum{iptc}(rs);
  cols(ncounter+(1:numel (rs))) = gnum{iptc}(cs);
  vals(ncounter+(1:numel (rs))) = vs;
  ncounter = ncounter + numel (rs);

  rhs_loc = op_f_v_tp (sp{iptc}, msh{iptc}, f);
  rhs(gnum{iptc}) = rhs(gnum{iptc}) + rhs_loc;
end

mat = sparse (rows, cols, vals, ndof, ndof);

% Apply Neumann boundary conditions
for iref = nmnn_sides
  for bnd_side = 1:boundaries(iref).nsides
    iptc = boundaries(iref).patches(bnd_side);
    iside = boundaries(iref).faces(bnd_side);
    msh_side = msh_eval_boundary_side (msh{iptc}, iside);
    sp_side  = sp_eval_boundary_side (sp, msh_side);

    x = squeeze (msh_side.geo_map(1,:,:));
    y = squeeze (msh_side.geo_map(2,:,:));
    gval = reshape (g (x, y, iref), 2, msh_side.nqn, msh_side.nel);
    rhs_nmnn = op_f_v (sp_side, msh_side, gval);
    global_dofs = gnum{iptc}(sp_side.dofs);
    rhs(global_dofs) = rhs(global_dofs) + rhs_nmnn;
  end
end

% Apply Dirichlet boundary conditions
u = zeros (ndof, 1);
[u_drchlt, drchlt_dofs] = mp_sp_drchlt_l2_proj (sp, msh, h, gnum, boundaries, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - mat(int_dofs, drchlt_dofs) * u_drchlt;

% Solve the linear system
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

end

%!demo
%! ex_plane_strain_Lshaped_mp
