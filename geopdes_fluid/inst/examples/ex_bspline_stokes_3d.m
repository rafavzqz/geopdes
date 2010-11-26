% EX_BSPLINE_STOKES_3D: Solve a Stokes flow problem on a three-dimensional domain.
%
% In order to make it work, a script file with the data must also be executed.
% The available data files for this problem are:
%
% - test_stokes_3d_symdrivcav
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

% Compute the geometry structure
geometry     = geo_load (geo_name);

% Compute the msh structure
breaks = {(linspace (0, 1, nbreaks(1))), (linspace (0, 1, nbreaks(2))), (linspace (0, 1, nbreaks(3)))};
knotsp = kntbrkdegreg (breaks, degree, regularity);
[qn, qw] = msh_set_quad_nodes (breaks, msh_gauss_nodes (nquad));
msh = msh_3d_tensor_product (breaks, qn, qw); 
msh = msh_push_forward_3d (msh, geometry);

% Compute the space structure
[spv, spp] = sp_bspline_th_3d_phys (knotsp, degree, msh); 

% Assemble the matrices
[x, y, z] = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)), squeeze (msh.geo_map(3,:,:)));
A = op_gradu_gradv (spv, spv, msh, mu (x, y, z)); 
B = op_div_v_q (spv, spp, msh); 
M = op_u_v (spp, spp, msh, ones (size (x))); 
E = sum (M, 1) / sum (sum (M)); 
F = op_f_v (spv, msh, f (x, y, z)); 

vel   = zeros (spv.ndof, 1);
press = zeros (spp.ndof, 1);

% Apply Dirichlet boundary conditions
[vel_drchlt, drchlt_dofs] = sp_drchlt_l2_proj(spv, msh, h, drchlt_sides);
vel(drchlt_dofs) = vel_drchlt;
int_dofs = setdiff (1:spv.ndof, drchlt_dofs);
nintdofs = numel (int_dofs);

% Solve the linear system
mat = [ A(int_dofs, int_dofs), -B(:,int_dofs).',            sparse(nintdofs, 1);
       -B(:,int_dofs),          sparse(spp.ndof, spp.ndof), E';
        sparse(1, nintdofs),    E,                          0];

rhs = [F(int_dofs)-A(int_dofs, drchlt_dofs)*vel(drchlt_dofs); 
       zeros(spp.ndof, 1); 
       0];

sol = mat \ rhs;

vel(int_dofs) = sol(1:nintdofs);
press = sol(1+nintdofs:end-1);

% Postprocessing
fprintf ('results being saved in: %s_velocity.vts and %s_pressure.vts\n', output_file, output_file)
sp_to_vtk_3d (vel, spv, geometry, vtk_pts, sprintf ('%s_velocity.vts', output_file), 'velocity')
sp_to_vtk_3d (press, spp, geometry, vtk_pts, sprintf ('%s_pressure.vts', output_file), 'pressure')

if (exist ('uex', 'var'))
  l2_error = sp_l2_error (spv, msh, vel, uex)
  if (exist ('graduex', 'var'))
    h1_error = sp_h1_error (spv, msh, vel, uex, graduex)
  end
end

%!demo
%! test_stokes_3d_symdrivcav
%! ex_bspline_stokes_3d
