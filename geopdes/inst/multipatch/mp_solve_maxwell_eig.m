% MP_SOLVE_MAXWELL_EIG: Solve the Maxwell eigenvalue problem with a B-spline discretization in a multipatch domain.
%
% Example to solve the problem
%
%    curl (1/mu(x) curl (u)) = lambda (epsilon(x) u)   in Omega
%      (1/mu(x) curl(u)) x n = 0                       on Gamma_N
%                      u x n = 0                       on Gamma_D
%
% where the domain \Omega is formed by several patches of the form F((0,1)^n).
%
% USAGE:
%
%  [geometry, msh, space, eigv, eigf, gnum, dofs_ornt] = 
%                  mp_solve_maxwell_eig (problem_data, method_data)
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
%  msh:       cell array of mesh objects (see msh_cartesian)
%  space:     cell array of space objects (see sp_vector_curl_transform)
%  eigv:      the computed eigenvalues
%  eigf:      degrees of freedom of the associated eigenfunctions
%  gnum:      global numbering of the degrees of freedom, for postprocessing
%               The sign indicates the global orientation of the vector.
%  dofs_ornt: orientation of the degrees of freedom (same as the sign of gnum)
%
% See also EX_MAXWELL_EIG_LSHAPED_MP for an example
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

function [geometry, msh, sp, eigv, eigf, gnum, dofs_ornt] = ...
              mp_solve_maxwell_eig (problem_data, method_data)

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
sp  = cell (1, npatch);
for iptc = 1:npatch
  [knots, zeta] = ...
         kntrefine (geometry(iptc).nurbs.knots, nsub-1, degree, regularity);
  [knots_hcurl, degree_hcurl] = knt_derham (knots, degree, 'Hcurl');

% Construct msh structure
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (zeta, rule);
  msh{iptc} = msh_cartesian (zeta, qn, qw, geometry(iptc));

% Construct space structure
  scalar_spaces = cell (msh{iptc}.ndim, 1);
  for idim = 1:msh{iptc}.ndim
    scalar_spaces{idim} = sp_bspline (knots_hcurl{idim}, degree_hcurl{idim}, msh{iptc});
  end
  sp{iptc} = sp_vector_curl_transform (scalar_spaces, msh{iptc});
  clear scalar_spaces
end

[gnum, ndof, dofs_ornt] = mp_interface_hcurl (interfaces, sp);

nc_stiff = 0; nc_mass = 0;
for iptc = 1:npatch

% Assemble the matrices setting the orientation
  if (msh{iptc}.rdim == 2)
    invmu = @(x,y) 1./c_magn_perm (x,y);
  elseif (msh{iptc}.rdim == 3)
    invmu = @(x,y,z) 1./c_magn_perm (x,y,z);
  end
  [rs, cs, vs] = op_curlu_curlv_tp (sp{iptc}, sp{iptc}, msh{iptc}, invmu);
  rows_stiff(nc_stiff+(1:numel (rs))) = gnum{iptc}(rs);
  cols_stiff(nc_stiff+(1:numel (rs))) = gnum{iptc}(cs);
  vs = dofs_ornt{iptc}(rs)' .* vs .* dofs_ornt{iptc}(cs)';
  vals_stiff(nc_stiff+(1:numel (rs))) = vs;
  nc_stiff = nc_stiff + numel (rs);

  [rs, cs, vs] = op_u_v_tp (sp{iptc}, sp{iptc}, msh{iptc}, c_elec_perm);
  rows_mass(nc_mass+(1:numel (rs))) = gnum{iptc}(rs);
  cols_mass(nc_mass+(1:numel (rs))) = gnum{iptc}(cs);
  vs = dofs_ornt{iptc}(rs)' .* vs .* dofs_ornt{iptc}(cs)';
  vals_mass(nc_mass+(1:numel (rs))) = vs;
  nc_mass = nc_mass + numel (rs);
end

clear rs cs vs
stiff_mat = sparse (rows_stiff, cols_stiff, vals_stiff, ndof, ndof);
clear rows_stiff cols_stiff vals_stiff
mass_mat  = sparse (rows_mass, cols_mass, vals_mass, ndof, ndof);
clear rows_mass cols_mass vals_mass

% Apply homogeneous Dirichlet boundary conditions
drchlt_dofs = [];
for iref = drchlt_sides
  for bnd_side = 1:boundaries(iref).nsides
    iptc = boundaries(iref).patches(bnd_side);
    iside = boundaries(iref).faces(bnd_side);
    global_dofs = gnum{iptc}(sp{iptc}.boundary(iside).dofs);
    drchlt_dofs = union (drchlt_dofs, global_dofs);
  end
end
int_dofs = setdiff (1:ndof, drchlt_dofs);

% Solve the eigenvalue problem
eigf = zeros (ndof, numel(int_dofs));
[eigf(int_dofs, :), eigv] = ...
    eig (full (stiff_mat(int_dofs, int_dofs)), full (mass_mat(int_dofs, int_dofs)));
eigv = sort (diag (eigv));

for iptc = 1:npatch
  gnum{iptc} = gnum{iptc} .*dofs_ornt{iptc};
end

end

%!demo
%! ex_maxwell_eig_Lshaped_mp

%!demo
%! ex_maxwell_eig_cube_mp

%!demo
%! ex_maxwell_eig_fichera_mp

%!demo
%! ex_maxwell_eig_thick_L_mp
