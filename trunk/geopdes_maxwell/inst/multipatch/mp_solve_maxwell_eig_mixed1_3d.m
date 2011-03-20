% MP_SOLVE_MAXWELL_EIG_MIXED1_3D: Solve the 3d Maxwell eigenvalue problem with a mixed formulation, and a B-spline discretization, in a multipatch domain.
%
% Example to solve the problem
%
%    curl (1/mu(x) curl (u)) = lambda (epsilon(x) u)   in Omega = F((0,1)^3)
%          div (epsilon(x) u) = 0                      in Omega
%       (1/mu(x) curl(u)) x n = 0                      on Gamma_N
%                       u x n = 0                      on Gamma_D
%
% with the variational mixed formulation
%
%    \int (1/mu(x) curl(u) curl(v)) + \int (epsilon(x) v grad(p))
%                = lambda \int (epsilon(x) u v),   \forall v \in H_0(curl),
%                                     \int (epsilon(x) u grad(q)) = 0,
%                                                  \forall q \in H^1_0.
%
% where the domain \Omega is formed by several patches of the form F((0,1)^3).
%
% USAGE:
%
%  [geometry, msh, space, sp_mul, eigv, eigf, gnum, dofs_ornt, gnum_mul] = 
%                  mp_solve_maxwell_eig_mixed1_3d (problem_data, method_data)
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
%    - n_sub:      number of subdivisions for refinement.
%    - nquad:      number of points for Gaussian quadrature rule
%
% OUTPUT:
%
%  geometry:  array of geometry structures (see mp_geo_load)
%  msh:       array of mesh structures (see msh_push_forward_3d)
%  space:     array of space structures (see sp_bspline_curl_transform_3d)
%  sp_mul:    array of space structures for the multiplier (see sp_bspline_3d_phys)
%  eigv:      the computed eigenvalues
%  eigf:      degrees of freedom of the associated eigenfunctions
%  gnum:      global numbering of the degrees of freedom, for postprocessing
%  dofs_ornt: orientation of the degrees of freedom, for postprocessing
%  gnum_mul   global numbering of the degrees of freedom of the multiplier
%
% See also EX_MAXWELL_EIG_MIXED1_CUBE_MP for an example
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

function [geometry, msh, sp, sp_mul, eigv, eigf, gnum, dofs_ornt, gnum_mul] = ...
              mp_solve_maxwell_eig_mixed1_3d (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} sprintf('= problem_data.(data_names{iopt});')]);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} sprintf('= method_data.(data_names{iopt});')]);
end

% Construct geometry structure
[geometry, boundaries, interfaces] = mp_geo_load (geo_name);
npatch = numel (geometry);

msh = cell (1, npatch); 
sp = cell (1, npatch);
for iptc = 1:npatch
  [knots, zeta] = ...
         kntrefine (geometry(iptc).nurbs.knots, n_sub, degree, regularity);
  [knots_u1, knots_u2, knots_u3, degree1, degree2, degree3] = ...
                                                 knt_derham (knots, degree);

% Construct msh structure
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (zeta, rule);
  msh{iptc} = msh_3d_tensor_product (zeta, qn, qw);
  msh{iptc} = msh_push_forward_3d (msh{iptc}, geometry(iptc));

% Construct space structure
  sp{iptc}     = sp_bspline_curl_transform_3d (knots_u1, knots_u2, knots_u3, ...
                                         degree1, degree2, degree3, msh{iptc});
  sp_mul{iptc} = sp_bspline_3d_phys (knots, degree, msh{iptc});
end

[gnum, ndof, dofs_ornt] = mp_interface_hcurl_3d (interfaces, sp);
[gnum_mul, ndof_mul] = mp_interface_3d (interfaces, sp_mul);

stiff_mat  = spalloc (ndof, ndof, ndof);
mass_mat   = spalloc (ndof, ndof, ndof);
saddle_mat = spalloc (ndof_mul, ndof, ndof_mul);

for iptc = 1:npatch

% Precompute the coefficients
  x = squeeze (msh{iptc}.geo_map(1,:,:));
  y = squeeze (msh{iptc}.geo_map(2,:,:));
  z = squeeze (msh{iptc}.geo_map(3,:,:));

  epsilon = reshape (c_elec_perm (x, y, z), msh{iptc}.nqn, msh{iptc}.nel);
  mu      = reshape (c_magn_perm (x, y, z), msh{iptc}.nqn, msh{iptc}.nel);

% Assemble the matrices
  stiff_loc  = op_curlu_curlv_3d (sp{iptc}, sp{iptc}, msh{iptc}, 1./mu);
  mass_loc   = op_u_v (sp{iptc}, sp{iptc}, msh{iptc}, epsilon);
  saddle_loc = op_v_gradp (sp{iptc}, sp_mul{iptc}, msh{iptc}, epsilon);

% Set the orientation
  ornt_matrix = spdiags (dofs_ornt{iptc}', 0, sp{iptc}.ndof, sp{iptc}.ndof);
  stiff_loc = ornt_matrix * stiff_loc * ornt_matrix;
  mass_loc = ornt_matrix * mass_loc * ornt_matrix;
  saddle_loc = saddle_loc * ornt_matrix;

  stiff_mat(gnum{iptc},gnum{iptc}) = stiff_mat(gnum{iptc},gnum{iptc}) + stiff_loc;
  mass_mat(gnum{iptc},gnum{iptc}) = mass_mat(gnum{iptc},gnum{iptc}) + mass_loc;
  saddle_mat(gnum_mul{iptc},gnum{iptc}) = ...
           saddle_mat(gnum_mul{iptc},gnum{iptc}) + saddle_loc;
end

% Apply homogeneous Dirichlet boundary conditions
drchlt_dofs = [];
drchlt_dofs_mul = [];
for iref = drchlt_sides
  for bnd_side = 1:boundaries(iref).nsides
    iptc = boundaries(iref).patches(bnd_side);
    iside = boundaries(iref).faces(bnd_side);
    global_dofs = gnum{iptc}(sp{iptc}.boundary(iside).dofs);
    drchlt_dofs = [drchlt_dofs global_dofs];
    global_dofs = gnum_mul{iptc}(sp_mul{iptc}.boundary(iside).dofs);
    drchlt_dofs_mul = [drchlt_dofs_mul global_dofs];
  end
end
drchlt_dofs = unique (drchlt_dofs);
drchlt_dofs_mul = unique (drchlt_dofs_mul);
int_dofs = setdiff (1:ndof, drchlt_dofs);
int_dofs_mul = setdiff (1:ndof_mul, drchlt_dofs_mul);

% Solve the eigenvalue problem
stiff_mat  = stiff_mat (int_dofs, int_dofs);
mass_mat   = mass_mat (int_dofs, int_dofs);
saddle_mat = saddle_mat (int_dofs_mul, int_dofs);

A = [stiff_mat, saddle_mat.'; ...
     saddle_mat, sparse(numel(int_dofs_mul),numel(int_dofs_mul))];
M = [mass_mat, sparse(numel(int_dofs), numel(int_dofs_mul)); ...
     sparse(numel(int_dofs_mul), numel(int_dofs)+numel(int_dofs_mul))];

[eigf, eigv] = eig (full(A), full(M));
eigv = diag (eigv);

end
