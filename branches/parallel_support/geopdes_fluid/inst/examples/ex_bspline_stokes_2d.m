% EX_BSPLINE_STOKES_2D: Solve a Stokes flow problem on a two-dimensional domain.
%
% In order to make it work, a script file with the data must also be executed.
% The available data files for this problem are:
%
% - test_stokes_square
% - test_stokes_square_bc
% - test_stokes_annulus
% - test_stokes_symdrivcav
%
% Copyright (C) 2009, 2010 Carlo de Falco
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
geometry = geo_load (geo_name);

% Construct msh structure
breaks = {(linspace (0, 1, nbreaks(1))), (linspace (0, 1, nbreaks(2)))};
knotsp = kntbrkdegreg (breaks, degree, regularity);
[qn, qw] = msh_set_quad_nodes (breaks, msh_gauss_nodes (nquad));
msh = msh_2d_tensor_product (breaks, qn, qw); 
msh = msh_push_forward_2d (msh, geometry, 'der2', der2);

% Construct space structure
[spv, spp] = feval (fun_space, knotsp, degree, msh);

% Assemble the matrices
[x, y] = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)));
A = op_gradu_gradv (spv, spv, msh, mu (x, y)); 
B = op_div_v_q (spv, spp, msh); 
M = op_u_v (spp, spp, msh, ones (size (x))); 
E = sum (M, 1) / sum (sum (M));
F = op_f_v (spv, msh, f (x, y));

vel   = zeros (spv.ndof, 1);
press = zeros (spp.ndof, 1);

% Apply Dirichlet boundary conditions
[vel_drchlt, drchlt_dofs] = sp_drchlt_l2_proj(spv, msh, h, drchlt_sides);
vel(drchlt_dofs) = vel_drchlt;
int_dofs = setdiff (1:spv.ndof, drchlt_dofs);
nintdofs = numel (int_dofs);

% Solve the linear system
mat = [A(int_dofs, int_dofs), -B(:,int_dofs).', sparse(nintdofs, 1);
       -B(:,int_dofs), sparse(spp.ndof, spp.ndof), E';
       sparse(1, nintdofs), E, 0];
rhs = [F(int_dofs)-A(int_dofs, drchlt_dofs)*vel(drchlt_dofs); 
       B(:, drchlt_dofs)*vel(drchlt_dofs); 
       0];

sol = mat \ rhs;

vel(int_dofs) = sol(1:nintdofs);
press = sol(1+nintdofs:end-1);

% Postprocessing
[eu, F] = sp_eval_2d (vel, spv, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

figure()
if (exist ('uex', 'var'))
  l2_error = sp_l2_error (spv, msh, vel, uex)
  if (exist ('graduex', 'var'))
    h1_error = sp_h1_error (spv, msh, vel, uex, graduex)
  end

  if (exist ('pex', 'var'))
    p_l2_error = sp_l2_error (spp, msh, press, pex)
  end

  subplot(1,2,2)
  euex = uex (X, Y);
  quiver (X, Y, squeeze(euex(1,:,:)), squeeze(euex(2,:,:)))
  axis equal
  title('Exact solution')
  
  subplot(1,2,1)
end

quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
axis equal
title('Computed solution')

fprintf ('results being saved in: %s_velocity.vts and %s_pressure.vts\n', output_file, output_file)
sp_to_vtk_2d (vel, spv, geometry, vtk_pts, sprintf ('%s_velocity.vts', output_file), 'velocity')
sp_to_vtk_2d (press, spp, geometry, vtk_pts, sprintf ('%s_pressure.vts', output_file), 'pressure')

%!demo
%! test_stokes_square
%! ex_bspline_stokes_2d

%!demo
%! test_stokes_square_bc
%! ex_bspline_stokes_2d

%!demo
%! test_stokes_symdrivcav
%! ex_bspline_stokes_2d

%!demo
%! test_stokes_annulus
%! ex_bspline_stokes_2d
