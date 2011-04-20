% EX_BSPLINE_STOKES_2D_MP: Solve a Stokes flow problem on a two-dimensional multipatch domain.
%
% Copyright (C) 2011 Carlo de Falco
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
npatch = numel (geometry);

ndofp = 0;
for iptc = 1:npatch
% Construct msh structure
  breaks = {(linspace (0, 1, nbreaks(1))), (linspace (0, 1, nbreaks(2)))};
  knotsp = kntbrkdegreg (breaks, degree, regularity);
  [qn, qw] = msh_set_quad_nodes (breaks, msh_gauss_nodes (nquad));
  msh{iptc} = msh_2d_tensor_product (breaks, qn, qw); 
  msh{iptc} = msh_push_forward_2d (msh{iptc}, geometry(iptc));

% Construct space structure
  [spv{iptc}, spp{iptc}] = sp_bspline_th_2d_phys (knotsp, degree, msh{iptc});

  gnump{iptc} = ndofp + (1:spp{iptc}.ndof);
  ndofp = ndofp + spp{iptc}.ndof;
end

% Create a correspondence between patches on the interfaces
[gnum, ndof] = mp_interface_vector_2d (interfaces, spv);

% Compute and assemble the matrices
A = spalloc (ndof, ndof, ndof);
M = spalloc (ndofp, ndofp, ndofp);
B = spalloc (ndofp, ndof, ndofp);
F = zeros (ndof, 1);

for iptc = 1:npatch
  [x, y] = deal (squeeze (msh{iptc}.geo_map(1,:,:)), squeeze (msh{iptc}.geo_map(2,:,:)));
  A_loc = op_gradu_gradv (spv{iptc}, spv{iptc}, msh{iptc}, mu (x, y));
  B_loc = op_div_v_q (spv{iptc}, spp{iptc}, msh{iptc});
  M_loc = op_u_v (spp{iptc}, spp{iptc}, msh{iptc}, ones (size (x))); 
  F_loc = op_f_v (spv{iptc}, msh{iptc}, f (x, y));

  A(gnum{iptc},gnum{iptc})   = A(gnum{iptc},gnum{iptc}) + A_loc;
  M(gnump{iptc},gnump{iptc}) = M(gnump{iptc},gnump{iptc}) + M_loc;
  B(gnump{iptc},gnum{iptc})  = B(gnump{iptc},gnum{iptc}) + B_loc;
  F(gnum{iptc}) = F(gnum{iptc}) + F_loc;
end
E = sum (M, 1) / sum (sum (M));

vel   = zeros (ndof, 1);
press = zeros (ndofp, 1);

% Apply Dirichlet boundary conditions
[vel_drchlt, drchlt_dofs] = mp_sp_drchlt_l2_proj (spv, msh, h, gnum, boundaries, drchlt_sides);
vel(drchlt_dofs) = vel_drchlt;

int_dofs = setdiff (1:ndof, drchlt_dofs);
nintdofs = numel (int_dofs);

% Solve the linear system
mat = [A(int_dofs, int_dofs), -B(:,int_dofs).', sparse(nintdofs, 1);
       -B(:,int_dofs), sparse(ndofp, ndofp), E';
       sparse(1, nintdofs), E, 0];
rhs = [F(int_dofs)-A(int_dofs, drchlt_dofs)*vel(drchlt_dofs); 
       B(:, drchlt_dofs)*vel(drchlt_dofs); 
       0];

sol = mat \ rhs;

vel(int_dofs) = sol(1:nintdofs);
press = sol(1+nintdofs:end-1);

% Postprocessing
if (exist ('uex', 'var'))
  for iptc = 1:npatch
    l2_error(iptc) = sp_l2_error (spv{iptc}, msh{iptc}, vel(gnum{iptc}), uex);
  end
  l2_error = sqrt (sum (l2_error .* l2_error))

  if (exist ('graduex', 'var'))
    for iptc = 1:npatch
      h1_error(iptc) = sp_h1_error (spv{iptc}, msh{iptc}, vel(gnum{iptc}), uex, graduex);
    end
    h1_error = sqrt (sum (h1_error .* h1_error))
  end

  if (exist ('pex', 'var'))
    for iptc = 1:npatch
      p_l2_error(iptc) = sp_l2_error (spp{iptc}, msh{iptc}, press(gnump{iptc}), pex);
    end
    p_l2_error = sqrt (sum (p_l2_error .* p_l2_error))
  end
end

fprintf ('results being saved in: %s_velocity.pvd and %s_pressure.pvd\n', output_file, output_file)
mp_sp_to_vtk_2d (vel, spv, geometry, gnum, vtk_pts, sprintf ('%s_velocity', output_file), 'velocity')
mp_sp_to_vtk_2d (press, spp, geometry, gnump, vtk_pts, sprintf ('%s_pressure', output_file), 'pressure')

