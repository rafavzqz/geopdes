% MP_SOLVE_MAXWELL_EIG_3D: Solve the 3d Maxwell eigenvalue problem with a B-spline discretization, in a multipatch domain.
%
% Example to solve the problem
%
%    curl (1/mu(x) curl (u)) = lambda (epsilon(x) u)   in Omega
%      (1/mu(x) curl(u)) x n = 0                       on Gamma_N
%                      u x n = 0                       on Gamma_D
%
% where the domain \Omega is formed by several patches of the form F((0,1)^3).
%
% USAGE:
%
%  [geometry, msh, space, eigv, eigf, gnum, dofs_ornt] = 
%                  mp_solve_maxwell_eig_3d (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_elec_perm:  electric permittivity (epsilon in the equation)
%    - c_magn_perm:  magnetic permeability (mu in the equation)
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
%  geometry:  array of geometry structures (see mp_geo_load)
%  msh:       array of mesh structures (see msh_push_forward_3d)
%  space:     array of space structures (see sp_bspline_curl_transform_3d)
%  eigv:      the computed eigenvalues
%  eigf:      degrees of freedom of the associated eigenfunctions
%  gnum:      global numbering of the degrees of freedom, for postprocessing
%  dofs_ornt: orientation of the degrees of freedom, for postprocessing
%
% See also EX_MAXWELL_EIG_CUBE_MP for an example
%
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

function [geometry, msh, sp, eigv, eigf, gnum, dofs_ornt] = ...
              mp_solve_maxwell_eig_3d (problem_data, method_data)

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
[geometry, boundaries, interfaces] = mp_geo_load (geo_name);
npatch = numel (geometry);

msh = cell (1, npatch); 
sp = cell (1, npatch);
for iptc = 1:npatch
  [knots, zeta] = ...
         kntrefine (geometry(iptc).nurbs.knots, nsub-1, degree, regularity);
  [knots_u1, knots_u2, knots_u3, degree1, degree2, degree3] = ...
                                                 knt_derham (knots, degree);

% Construct msh structure
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (zeta, rule);
  msh{iptc} = msh_3d_tensor_product (zeta, qn, qw);
  msh{iptc} = msh_push_forward_3d (msh{iptc}, geometry(iptc));

% Construct space structure
  sp{iptc} = sp_bspline_curl_transform_3d (knots_u1, knots_u2, knots_u3, ...
                                       degree1, degree2, degree3, msh{iptc});
end

[gnum, ndof, dofs_ornt] = mp_interface_hcurl_3d (interfaces, sp);

nent = sum (cellfun (@(x, y, z) x.nel * y.nsh_max * z.nsh_max, msh, sp, sp));
rows = zeros (nent, 1);
cols = zeros (nent, 1);
vals_stiff = zeros (nent, 1);
vals_mass  = zeros (nent, 1);

ncounter = 0;
for iptc = 1:npatch

% Precompute the coefficients
  x = squeeze (msh{iptc}.geo_map(1,:,:));
  y = squeeze (msh{iptc}.geo_map(2,:,:));
  z = squeeze (msh{iptc}.geo_map(3,:,:));

  epsilon = reshape (c_elec_perm (x, y, z), msh{iptc}.nqn, msh{iptc}.nel);
  mu      = reshape (c_magn_perm (x, y, z), msh{iptc}.nqn, msh{iptc}.nel);

% Assemble the matrices setting the orientation
  [rs, cs, vs] = op_curlu_curlv_3d (sp{iptc}, sp{iptc}, msh{iptc}, 1./mu);
  rows(ncounter+(1:numel (rs))) = gnum{iptc}(rs);
  cols(ncounter+(1:numel (rs))) = gnum{iptc}(cs);
  vs = dofs_ornt{iptc}(rs)' .* vs .* dofs_ornt{iptc}(cs)';
  vals_stiff(ncounter+(1:numel (rs))) = vs;

  [rs, cs, vs] = op_u_v (sp{iptc}, sp{iptc}, msh{iptc}, epsilon);
  vs = dofs_ornt{iptc}(rs)' .* vs .* dofs_ornt{iptc}(cs)';
  vals_mass(ncounter+(1:numel (rs))) = vs;
  ncounter = ncounter + numel (rs);
end

stiff_mat = sparse (rows, cols, vals_stiff, ndof, ndof);
mass_mat  = sparse (rows, cols, vals_mass, ndof, ndof);

% Apply homogeneous Dirichlet boundary conditions
drchlt_dofs = [];
for iref = drchlt_sides
  for bnd_side = 1:boundaries(iref).nsides
    iptc = boundaries(iref).patches(bnd_side);
    iside = boundaries(iref).faces(bnd_side);
    global_dofs = gnum{iptc}(sp{iptc}.boundary(iside).dofs);
    drchlt_dofs = [drchlt_dofs global_dofs];
  end
end
drchlt_dofs = unique (drchlt_dofs);
int_dofs = setdiff (1:ndof, drchlt_dofs);

% Solve the eigenvalue problem
eigf = zeros (ndof, numel(int_dofs));
[eigf(int_dofs, :), eigv] = ...
    eig (full (stiff_mat(int_dofs, int_dofs)), full (mass_mat(int_dofs, int_dofs)));
eigv = diag (eigv);

end

%!demo
%! ex_maxwell_eig_cube_mp

%!demo
%! ex_maxwell_eig_thick_L_mp

%!demo
%! ex_maxwell_eig_fichera_mp
