% EX_NURBS_LAPLACE_2D: Solve a 2d Laplace problem with a NURBS discretization. 
%
% Example to solve the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega = F((0,1)^2)
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% In order to make it work, a script file with the data must also be executed.
% The available data files for this problem are:
%
% - test_square
% - test_ring
% - test_ring_mixed_bc
% - test_plate_mixed_bc
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
degelev  = max (degree - (geometry.nurbs.order-1), 0);
nurbs    = nrbdegelev (geometry.nurbs, degelev);
[rknots, zeta, nknots] = kntrefine (nurbs.knots, n_sub, nurbs.order-1, regularity);

nurbs = nrbkntins (nurbs, nknots);
geometry = geo_load (nurbs);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (geometry.nurbs.knots, rule);
msh      = msh_2d_tensor_product (geometry.nurbs.knots, qn, qw);
msh      = msh_push_forward_2d (msh, geometry);
  
% Construct space structure
sp       = sp_nurbs_2d_phys (geometry.nurbs, msh);

% Precompute the coefficients
x = squeeze (msh.geo_map(1,:,:));
y = squeeze (msh.geo_map(2,:,:));
  
epsilon = reshape (c_diff (x, y), msh.nqn, msh.nel);
fval    = reshape (f (x, y), msh.nqn, msh.nel) ;
 
% Assemble the matrices
stiff_mat = op_gradu_gradv (sp, sp, msh, epsilon);
rhs       = op_f_v (sp, msh, fval);

% Apply Neumann boundary conditions
for iside = nmnn_sides
  x = squeeze (msh.boundary(iside).geo_map(1,:,:));
  y = squeeze (msh.boundary(iside).geo_map(2,:,:));
  gval = reshape (g (x, y, iside), msh.boundary(iside).nqn, msh.boundary(iside).nel);

  rhs(sp.boundary(iside).dofs) = rhs(sp.boundary(iside).dofs) + ...
      op_f_v (sp.boundary(iside), msh.boundary(iside), gval);
end

% Apply Dirichlet boundary conditions
u = zeros (sp.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj(sp, msh, h, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:sp.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system
u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);

% Postprocessing
[eu, F] = sp_eval_2d (u, sp, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

if (exist ('uex', 'var'))
  error_l2 = sp_l2_error (sp, msh, u, uex)
  if (exist ('graduex', 'var'))
    error_h1 = sp_h1_error (sp, msh, u, uex, graduex)
  end

  subplot (1,2,2)
  surf (X, Y, uex (X,Y))
  title ('exact solution'), axis tight
  subplot(1,2,1)
end

surf (X, Y, eu);
title ('numerical solution'), axis tight

fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk_2d (u, sp, geometry, vtk_pts, output_file, 'u')


%!demo
%! test_ring_mixed_bc
%! ex_nurbs_laplace_2d

%!demo
%! test_square
%! ex_nurbs_laplace_2d

%!demo
%! test_ring
%! ex_nurbs_laplace_2d

%!demo
%! test_plate_mixed_bc
%! ex_nurbs_laplace_2d
