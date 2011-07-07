% MP_SOLVE_LAPLACE_3D: solve the Laplacian problem in a multipatch geometry.
%
% Example to solve the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% where the domain \Omega is formed by several patches of the form F((0,1)^3).
%
% USAGE:
%
%  [geometry, msh, space, u, gnum] = 
%          mp_solve_laplace_3d (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_diff:       diffusion coefficient (epsilon in the equation)
%    - f:            source term
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
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
%  msh:      cell array of mesh objects (see msh_3d)
%  space:    cell array of space objects (see sp_bspline_3d)
%  u:        the computed degrees of freedom
%  gnum:     global numbering of the degrees of freedom
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2010, 2011 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 2 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [geometry, msh, sp, u, gnum] = ...
              mp_solve_laplace_3d (problem_data, method_data)

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

msh = cell (1, npatch); 
sp = cell (1, npatch);
for iptc = 1:npatch

% Define the refined mesh, with tensor product structure
  [knots{iptc}, zeta{iptc}] = ...
         kntrefine (geometry(iptc).nurbs.knots, nsub-1, degree, regularity);

% Compute the quadrature rule
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (zeta{iptc}, rule);
  msh{iptc} = msh_3d (zeta{iptc}, qn, qw, geometry(iptc));

% Evaluate the discrete space basis functions in the quadrature points
  sp{iptc} = sp_bspline_3d (knots{iptc}, degree, msh{iptc});
end

% Create a correspondence between patches on the interfaces
[gnum, ndof] = mp_interface_3d (interfaces, sp);

% Compute and assemble the matrices 
rhs = zeros (ndof, 1);

ncounter = 0;
for iptc = 1:npatch
  [rs, cs, vs] = op_gradu_gradv_tp (sp{iptc}, sp{iptc}, msh{iptc}, c_diff);
  rows(ncounter+(1:numel (rs))) = gnum{iptc}(rs);
  cols(ncounter+(1:numel (rs))) = gnum{iptc}(cs);
  vals(ncounter+(1:numel (rs))) = vs;
  ncounter = ncounter + numel (rs);

  rhs_loc = op_f_v_tp (sp{iptc}, msh{iptc}, f);
  rhs(gnum{iptc}) = rhs(gnum{iptc}) + rhs_loc;
end

stiff_mat = sparse (rows, cols, vals, ndof, ndof);

for iref = nmnn_sides
  for bnd_side = 1:boundaries(iref).nsides
    iptc = boundaries(iref).patches(bnd_side);
    iside = boundaries(iref).faces(bnd_side);
    msh_side = msh_eval_boundary_side (msh{iptc}, iside);
    sp_side  = sp_eval_boundary_side (sp{iptc}, msh_side);

    x = squeeze (msh_side.geo_map(1,:,:));
    y = squeeze (msh_side.geo_map(2,:,:));
    z = squeeze (msh_side.geo_map(3,:,:));
    gval = reshape (g (x, y, z, iref), msh_side.nqn, msh_side.nel);
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
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system
u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);

end

%!demo
%! ex_laplace_cube_mp

%!demo
%! ex_laplace_thick_L_mp