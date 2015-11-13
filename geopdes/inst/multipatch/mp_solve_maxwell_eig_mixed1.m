% MP_SOLVE_MAXWELL_EIG_MIXED1: Solve the Maxwell eigenvalue problem with a mixed formulation, and a B-spline discretization, in a multipatch domain.
%
% Example to solve the problem
%
%    curl (1/mu(x) curl (u)) = lambda (epsilon(x) u)   in Omega 
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
% where the domain \Omega is formed by several patches of the form F((0,1)^n).
%
% USAGE:
%
%  [geometry, msh, space, sp_mul, eigv, eigf] = 
%                  mp_solve_maxwell_eig_mixed1 (problem_data, method_data)
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
%  geometry: array of geometry structures (see mp_geo_load)
%  msh:      multipatch mesh, consisting of several Cartesian meshes (see msh_multipatch)
%  space:    multipatch space, formed by several tensor product spaces plus the connectivity and orientation (see sp_multipatch)
%  sp_mul:   multipatch space for the multiplier (see sp_multipatch)
%  eigv:     the computed eigenvalues
%  eigf:     degrees of freedom of the associated eigenfunctions
%
% See also EX_MAXWELL_EIG_MIXED1_LSHAPED_MP for an example
%
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

function [geometry, msh, space, space_mul, eigv, eigf] = ...
              mp_solve_maxwell_eig_mixed1 (problem_data, method_data)

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
[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (geo_name);
npatch = numel (geometry);

msh = cell (1, npatch);
sp  = cell (1, npatch);
sp_mul = cell (1, npatch);
for iptc = 1:npatch
  [knots, zeta] = ...
         kntrefine (geometry(iptc).nurbs.knots, nsub-1, degree, regularity);
  [knots_hcurl, degree_hcurl] = knt_derham (knots, degree, 'Hcurl');

% Construct msh structure
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (zeta, rule);
  msh{iptc} = msh_cartesian (zeta, qn, qw, geometry(iptc));

% Construct space structures
  scalar_spaces = cell (msh{iptc}.ndim, 1);
  for idim = 1:msh{iptc}.ndim
    scalar_spaces{idim} = sp_bspline (knots_hcurl{idim}, degree_hcurl{idim}, msh{iptc});
  end
  sp{iptc} = sp_vector (scalar_spaces, msh{iptc}, 'curl-preserving');
  clear scalar_spaces
  sp_mul{iptc} = sp_bspline (knots, degree, msh{iptc});
end

msh = msh_multipatch (msh, boundaries);
space = sp_multipatch (sp, msh, interfaces, boundary_interfaces);
space_mul = sp_multipatch (sp_mul, msh, interfaces, boundary_interfaces);
clear sp sp_mul

% Assemble the matrices setting the orientation
if (msh.rdim == 2)
  invmu = @(x,y) 1./c_magn_perm (x,y);
elseif (msh.rdim == 3)
  invmu = @(x,y,z) 1./c_magn_perm (x,y,z);
end
stiff_mat = op_curlu_curlv_mp (space, space, msh, invmu);
mass_mat = op_u_v_mp (space, space, msh, c_elec_perm);
saddle_mat = op_v_gradp_mp (space, space_mul, msh, c_elec_perm);

% Apply homogeneous Dirichlet boundary conditions
Nbnd = cumsum ([0, boundaries.nsides]);
boundary_gnum = space.boundary.gnum;
bnd_dofs = [];
boundary_gnum_mul = space_mul.boundary.gnum;
bnd_dofs_mul = [];
for iref = drchlt_sides
  iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
  bnd_dofs = union (bnd_dofs, [boundary_gnum{iref_patch_list}]);
  bnd_dofs_mul = union (bnd_dofs_mul, [boundary_gnum_mul{iref_patch_list}]);
end
drchlt_dofs = space.boundary.dofs (bnd_dofs);
drchlt_dofs_mul = space_mul.boundary.dofs (bnd_dofs_mul);
if (isempty(drchlt_dofs_mul))
  drchlt_dofs_mul = ndof_mul;
end

int_dofs = setdiff (1:space.ndof, drchlt_dofs);
int_dofs_mul = setdiff (1:space_mul.ndof, drchlt_dofs_mul);

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

%!demo
%! ex_maxwell_eig_mixed1_Lshaped_mp
