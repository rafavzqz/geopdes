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

subdomains = [2 3 4 5];
[geometry, interfaces, boundaries] = mp_extract_subdomains (geometry, interfaces, boundaries, subdomains);
npatch = numel (geometry);

for iptc = 1:npatch
  degelev  = max (degree - (geometry(iptc).nurbs.order-1), 0);
  nurbs    = nrbdegelev (geometry(iptc).nurbs, degelev);
  [rknots, zeta, nknots] = kntrefine (nurbs.knots, n_sub, nurbs.order-1, regularity);

  nurbs    = nrbkntins (nurbs, nknots);
  geometry(iptc) = geo_load (nurbs);

% Construct msh structure
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (geometry(iptc).nurbs.knots, rule);
  msh{iptc} = msh_3d_tensor_product (geometry(iptc).nurbs.knots, qn, qw);
  msh{iptc} = msh_push_forward_3d (msh{iptc}, geometry(iptc));

% Construct space structure
  sp_scalar = sp_nurbs_3d_phys (geometry(iptc).nurbs, msh{iptc});
  sp{iptc} = sp_scalar_to_vector_3d (sp_scalar, sp_scalar, sp_scalar, ...
                                     msh{iptc}, 'divergence', true);
end

% Create a correspondence between patches on the interfaces
[gnum, ndof] = mp_interface_vector_3d (interfaces, sp);

% Compute and assemble the matrices
mat = spalloc (ndof, ndof, ndof);
rhs = zeros (ndof, 1);

for iptc = 1:npatch
  [x, y, z] = deal (squeeze (msh{iptc}.geo_map(1,:,:)), squeeze (msh{iptc}.geo_map(2,:,:)), squeeze (msh{iptc}.geo_map(3,:,:)));
  
  mat_loc = op_su_ev (sp{iptc}, sp{iptc}, msh{iptc}, lam (x, y, z), mu (x, y, z)); 
  rhs_loc = op_f_v (sp{iptc}, msh{iptc}, f (x, y, z));

  mat(gnum{iptc},gnum{iptc}) = mat(gnum{iptc},gnum{iptc}) + mat_loc;
  rhs(gnum{iptc}) = rhs(gnum{iptc}) + rhs_loc;
end

% Apply Neumann boundary conditions
for iref = nmnn_sides
  for bnd_side = 1:boundaries(iref).nsides
    iptc = boundaries(iref).patches(bnd_side);
    iside = boundaries(iref).faces(bnd_side);
    x = squeeze (msh.boundary(iside).geo_map(1,:,:));
    y = squeeze (msh.boundary(iside).geo_map(2,:,:));
    z = squeeze (msh.boundary(iside).geo_map(3,:,:));
    gval = reshape (g (x, y, z, iref), ...
           msh{iptc}.boundary(iside).nqn, msh{iptc}.boundary(iside).nel);
    rhs_nmnn = ...
           op_f_v (sp{iptc}.boundary(iside), msh{iptc}.boundary(iside), gval);
    global_dofs = gnum{iptc}(sp{iptc}.boundary(iside).dofs);
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

% Postprocessing
if (exist ('uex', 'var'))
  for iptc = 1:npatch
    error_l2(iptc) = sp_l2_error (sp{iptc}, msh{iptc}, u(gnum{iptc}), uex);
  end
  error_l2 = sqrt (sum (error_l2 .* error_l2))

  if (exist ('graduex', 'var'))
    for iptc = 1:npatch
      error_h1(iptc) = sp_h1_error (sp{iptc}, msh{iptc}, u(gnum{iptc}), uex, graduex);
    end
    error_h1 = sqrt (sum (error_h1 .* error_h1))
  end
end

fprintf ('results being saved in: %s_displacement.pvd\n', output_file)
mp_sp_to_vtk_3d (u, sp, geometry, gnum, vtk_pts, sprintf ('%s_displacement', output_file), 'displacement')

%!demo
%! test_waveguide_3d_mp
%! ex_nurbs_linear_elasticity_3d_extract_mp
