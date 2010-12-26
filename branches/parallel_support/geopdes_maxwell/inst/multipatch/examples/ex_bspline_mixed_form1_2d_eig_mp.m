% EX_BSPLINE_MIXED_FORM1_2D_EIG_MP: Solve the 2d Maxwell eigenvalue problem with a mixed formulation, and a B-spline discretization, in a multipatch domain.
%
% Example to solve the problem
%
%    curl ( 1/mu(x) curl (u)) = lambda (epsilon(x) u)   in Omega 
%          div (epsilon(x) u) = 0                       in Omega
%       (1/mu(x) curl(u)) x n = 0                       on Gamma_N
%                       u x n = 0                       on Gamma_D
%
% with the variational mixed formulation
%
%    \int (1/mu(x) curl(u) curl(v)) + \int (epsilon(x) v grad(p))
%                = lambda \int (epsilon(x) u v),   \forall v \in H_0(curl),
%                                     \int (epsilon(x) u grad(q)) = 0,
%                                                  \forall q \in H^1_0.
%
% where the domain \Omega is formed by several patches of the form F((0,1)^2).
%
% In order to make it work, a script file with the data must also be executed.
% The available data files for this problem are:
%
% - test_maxwell_Lshaped_eig_mp
%
% Copyright (C) 2010 Rafael Vazquez
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

% Construct geometry structure
[geometry, boundaries, interfaces] = mp_geo_load (geo_name);
npatch = numel (geometry);

msh = cell (1, npatch); 
sp = cell (1, npatch);
for iptc = 1:npatch
  [knots, zeta] = ...
         kntrefine (geometry(iptc).nurbs.knots, n_sub, degree, regularity);
  [knots_u1, knots_u2, degree1, degree2] = knt_derham (knots, degree);

% Construct msh structure
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (zeta, rule);
  msh{iptc} = msh_2d_tensor_product (zeta, qn, qw);
  msh{iptc} = msh_push_forward_2d (msh{iptc}, geometry(iptc));

% Construct space structure
  sp{iptc} = sp_bspline_curl_transform_2d (knots_u1, knots_u2, ...
                                           degree1, degree2, msh{iptc});
  sp_mul{iptc} = sp_bspline_2d_phys (knots, degree, msh{iptc});
end

[gnum, ndof, dofs_ornt] = mp_interface_hcurl_2d (interfaces, sp);
[gnum_mul, ndof_mul] = mp_interface_2d (interfaces, sp_mul);

stiff_mat  = spalloc (ndof, ndof, ndof);
mass_mat   = spalloc (ndof, ndof, ndof);
saddle_mat = spalloc (ndof_mul, ndof, ndof_mul);

for iptc = 1:npatch

% Precompute the coefficients
  x = squeeze (msh{iptc}.geo_map(1,:,:));
  y = squeeze (msh{iptc}.geo_map(2,:,:));

  epsilon = reshape (c_elec_perm (x, y), msh{iptc}.nqn, msh{iptc}.nel);
  mu      = reshape (c_magn_perm (x, y), msh{iptc}.nqn, msh{iptc}.nel);

% Assemble the matrices
  stiff_loc  = op_curlu_curlv_2d (sp{iptc}, sp{iptc}, msh{iptc}, 1./mu);
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

eigv = eig (full(A), full(M));
eigv = sort (eigv);


fprintf ('First computed eigenvalues: \n')
disp (eigv(1:5))

%!demo
%! test_maxwell_Lshaped_eig_mp
%! ex_bspline_mixed_form1_2d_eig_mp

