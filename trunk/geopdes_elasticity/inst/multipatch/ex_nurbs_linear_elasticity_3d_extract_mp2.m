% EX_NURBS_LINEAR_ELASTICITY_3D_EXTRACT_MP: Solve a linear elasticity problem on 
% a portion of a 3-dimensional multipatch domain.
%
% The problem solved is
%
%      - div (sigma(u)) = f    in Omega
%      sigma(u) \cdot n = g    on Gamma_N
%                     u = h    on Gamma_D
%
% with   sigma(u) = mu*(grad(u) + grad(u)^t) + lambda*div(u)*I,
% and the domain \Omega is formed by several patches of the form F((0,1)^3).
%
%   u:          displacement vector
%   sigma:      Cauchy stress tensor
%   lambda, mu: Lame' parameters
%   I:          identity tensor
%
% In order to make it work, a script file with the data must also be executed.
% The available data files for this problem are:
%
% - test_waveguide_3d_mp
%
% Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
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

% Construct geometry structure, and information for interfaces and boundaries
[geometry, boundaries, interfaces] = mp_geo_load (geo_name);

%solve problem for the external shell
subdomains = [2:5 7:10 12:15];
[geo_ext, int_ext, bnd_ext] = mp_extract_subdomains2 (geometry, interfaces, boundaries, subdomains);
npatch = numel (geo_ext);

for iptc = 1:npatch
  degelev  = max (degree - (geo_ext(iptc).nurbs.order-1), 0);
  nurbs    = nrbdegelev (geo_ext(iptc).nurbs, degelev);
  [rknots, zeta, nknots] = kntrefine (nurbs.knots, n_sub, nurbs.order-1, regularity);

  nurbs    = nrbkntins (nurbs, nknots);
  geo_ext(iptc) = geo_load (nurbs);

  % Construct msh structure
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (geo_ext(iptc).nurbs.knots, rule);
  msh_ext{iptc} = msh_3d_tensor_product (geo_ext(iptc).nurbs.knots, qn, qw);
  msh_ext{iptc} = msh_push_forward_3d (msh_ext{iptc}, geo_ext(iptc));

  % Construct space structure
  sp_scalar = sp_nurbs_3d_phys (geo_ext(iptc).nurbs, msh_ext{iptc});
  sp_ext{iptc} = sp_scalar_to_vector_3d (sp_scalar, sp_scalar, sp_scalar, ...
                                         msh_ext{iptc}, 'divergence', true);
end

% Create a correspondence between patches on the interfaces
[gnum_ext, ndof_ext, gnum_ext_bnd] = mp_interface_vector_3d_sub (int_ext, bnd_ext, sp_ext); 

% Compute and assemble the matrices
mat = spalloc (ndof_ext, ndof_ext, ndof_ext);
rhs = zeros (ndof_ext, 1);

for iptc = 1:npatch
  [x, y, z] = deal (squeeze (msh_ext{iptc}.geo_map(1,:,:)), squeeze (msh_ext{iptc}.geo_map(2,:,:)), squeeze (msh_ext{iptc}.geo_map(3,:,:)));
  
  mat_loc = op_su_ev (sp_ext{iptc}, sp_ext{iptc}, msh_ext{iptc}, lam (x, y, z), mu (x, y, z)); 
  rhs_loc = op_f_v (sp_ext{iptc}, msh_ext{iptc}, f (x, y, z));

  mat(gnum_ext{iptc},gnum_ext{iptc}) = mat(gnum_ext{iptc},gnum_ext{iptc}) + mat_loc;
  rhs(gnum_ext{iptc}) = rhs(gnum_ext{iptc}) + rhs_loc;
end

% Apply Dirichlet boundary conditions
u = zeros (ndof_ext, 1);
[u_drchlt, drchlt_dofs] = mp_sp_drchlt_l2_proj (sp_ext, msh_ext, h, gnum_ext, bnd_ext, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:ndof_ext, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - mat(int_dofs, drchlt_dofs) * u_drchlt;

% Solve the linear system
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

fprintf ('results being saved in: %s_displacement.pvd\n', output_file)
mp_sp_to_vtk_3d (u, sp_ext, geo_ext, gnum_ext, vtk_pts, sprintf ('%s_displacement', output_file), 'displacement')

%solve problem for the internal domain
subdomains = [1 6 11];
[geo_int, int_int, bnd_int] = mp_extract_subdomains2 (geometry, interfaces, boundaries, subdomains);
npatch = numel (geo_int);

for iptc = 1:npatch
  degelev  = max (degree - (geo_int(iptc).nurbs.order-1), 0);
  nurbs    = nrbdegelev (geo_int(iptc).nurbs, degelev);
  [rknots, zeta, nknots] = kntrefine (nurbs.knots, n_sub, nurbs.order-1, regularity);

  nurbs    = nrbkntins (nurbs, nknots);
  geo_int(iptc) = geo_load (nurbs);

  % Construct msh structure
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (geo_int(iptc).nurbs.knots, rule);
  msh_int{iptc} = msh_3d_tensor_product (geo_int(iptc).nurbs.knots, qn, qw);
  msh_int{iptc} = msh_push_forward_3d (msh_int{iptc}, geo_int(iptc));

  % Construct space structure
  sp_scalar = sp_nurbs_3d_phys (geo_int(iptc).nurbs, msh_int{iptc});
  sp_int{iptc} = sp_scalar_to_vector_3d (sp_scalar, sp_scalar, sp_scalar, ...
                                         msh_int{iptc}, 'divergence', true);
end

% Create a correspondence between patches on the interfaces
[gnum_int, ndof_int, gnum_int_bnd] = mp_interface_vector_3d_sub (int_int, bnd_int, sp_int); 

% Compute and assemble the matrices
mat = spalloc (ndof_int, ndof_int, ndof_int);
rhs = zeros (ndof_int, 1);

for iptc = 1:npatch
  [x, y, z] = deal (squeeze (msh_int{iptc}.geo_map(1,:,:)), squeeze (msh_int{iptc}.geo_map(2,:,:)), squeeze (msh_int{iptc}.geo_map(3,:,:)));
  
  mat_loc = op_su_ev (sp_int{iptc}, sp_int{iptc}, msh_int{iptc}, lam (x, y, z), mu (x, y, z)); 
  rhs_loc = op_f_v (sp_int{iptc}, msh_int{iptc}, f (x, y, z));

  mat(gnum_int{iptc},gnum_int{iptc}) = mat(gnum_int{iptc},gnum_int{iptc}) + mat_loc;
  rhs(gnum_int{iptc}) = rhs(gnum_int{iptc}) + rhs_loc;
end

% Apply Dirichlet boundary conditions
uint = zeros (ndof_int, 1);
drchlt_dofs = [];
for jjj = [3:14]
  uint(gnum_int_bnd{jjj}{1}) = u(gnum_ext_bnd{jjj}{1});
  drchlt_dofs = union (drchlt_dofs, gnum_int_bnd{jjj}{1});
end

int_dofs = setdiff (1:ndof_int, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - mat(int_dofs, drchlt_dofs) * uint(drchlt_dofs);

% Solve the linear system
uint(int_dofs) = 0;%mat(int_dofs, int_dofs) \ rhs(int_dofs);

fprintf ('results being saved in: %s_int_displacement.pvd\n', output_file)
mp_sp_to_vtk_3d (uint, sp_int, geo_int, gnum_int, vtk_pts, sprintf ('%s_int_displacement', output_file), 'displacement')



%!demo
%! test_waveguide_3d_mp
%! ex_nurbs_linear_elasticity_3d_extract_mp
