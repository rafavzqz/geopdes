% EX_BSPLINE_LAPLACE_3D_MP: solve the Laplacian problem in a multipatch geometry.
%
% Example to solve the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% where the domain \Omega is formed by several patches of the form F((0,1)^3).
%
% In order to make it work, a script file with the data must also be executed.
% The available data files for this problem are:
%
% - test_cube_mp
% - test_thick_Lshaped_mp
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2010 Rafael Vazquez
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

% Construct geometry structure, and information for interfaces and boundaries
[geometry, boundaries, interfaces] = mp_geo_load (geo_name);
npatch = numel (geometry);

msh = cell (1, npatch); 
sp = cell (1, npatch);
for iptc = 1:npatch

% Define the refined mesh, with tensor product structure
  [knots{iptc}, zeta{iptc}] = ...
         kntrefine (geometry(iptc).nurbs.knots, n_sub, degree, regularity);

% Compute the quadrature rule
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (zeta{iptc}, rule);
  msh{iptc} = msh_3d_tensor_product (zeta{iptc}, qn, qw);
  msh{iptc} = msh_push_forward_3d (msh{iptc}, geometry(iptc));

% Evaluate the discrete space basis functions in the quadrature points
  sp{iptc} = sp_bspline_3d_phys (knots{iptc}, degree, msh{iptc});
end

% Create a correspondence between patches on the interfaces
[gnum, ndof] = mp_interface_3d (interfaces, sp);

% Compute and assemble the matrices 
stiff_mat = spalloc (ndof, ndof, ndof);
rhs = zeros (ndof, 1);
for iptc = 1:npatch

  x = squeeze (msh{iptc}.geo_map (1,:,:));
  y = squeeze (msh{iptc}.geo_map (2,:,:));
  z = squeeze (msh{iptc}.geo_map (3,:,:));

  epsilon = reshape (c_diff (x, y, z), msh{iptc}.nqn, msh{iptc}.nel);
  fval    = reshape (f(x, y, z), msh{iptc}.nqn, msh{iptc}.nel);

  A_loc   = op_gradu_gradv (sp{iptc}, sp{iptc}, msh{iptc}, epsilon);
  rhs_loc = op_f_v (sp{iptc}, msh{iptc}, fval);

  stiff_mat(gnum{iptc},gnum{iptc}) = stiff_mat(gnum{iptc},gnum{iptc}) + A_loc;
  rhs(gnum{iptc}) = rhs(gnum{iptc}) + rhs_loc;
end

for iref = nmnn_sides
  for bnd_side = 1:boundaries(iref).nsides
    iptc = boundaries(iref).patches(bnd_side);
    iside = boundaries(iref).faces(bnd_side);
    x = squeeze (msh{iptc}.boundary(iside).geo_map(1,:,:));
    y = squeeze (msh{iptc}.boundary(iside).geo_map(2,:,:));
    z = squeeze (msh{iptc}.boundary(iside).geo_map(3,:,:));
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
[u_drchlt, drchlt_dofs] = mp_sp_drchlt_l2_proj(sp, msh, h, gnum, boundaries, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system
u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);

% Compute the errors
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

% Postprocessing
fprintf ('The result is saved in the file %s.pvd \n \n', output_file);
mp_sp_to_vtk_3d (u, sp, geometry, gnum, vtk_pts, output_file, 'u')

%!demo
%! test_thick_Lshaped_mp
%! ex_bspline_laplace_3d_mp

%!demo
%! test_cube_mp
%! ex_bspline_laplace_3d_mp
