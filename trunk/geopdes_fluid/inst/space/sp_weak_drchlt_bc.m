% SP_WEAK_DRCHLT_BC: compute the matrix and right hand-side to impose
% the Dirichlet boundary conditions in weak form for the tangetial
% component, and in strong form for the normal component. To be used with
% the 'RT' spaces.
%
% The code computes the following terms in the left hand-side
% 
%  - \int_{Gamma_D} mu*{(\grad u)n \cdot v - (\grad v)n \cdot n + (Cpen /  he) * (u cdot v)}
%
% and in the right hand-side
%
%  - \int_{Gamma_D} mu*{(\grad v)n \cdot g + (Cpen / he) * (v cdot g)}
%
% with u the trial function, v the test function, he the normal characteristic length, 
%  and g the boundary condition to be imposed.
%
%
%   [N_mat, N_rhs, vel, normal_dofs] = solve_nitsche_weak_bc (space, msh, geometry, der2, bnd_sides, bnd_func, coeff, Cpen)
%
% INPUTS:
%     
%    space:     space object (see sp_vector_piola_transform_2d)
%    msh:       mesh object (see msh_2d)
%    geometry:  geometry object (see geo_load)
%    der2:      a logical telling whether to compute the second derivative.
%               For Raviart-Thomas (div-conforming) elements, it should be true.
%    bnd_sides: boundary sides on which the Dirichlet condition is imposed
%    bnd_func:  the condition to be imposed (g in the equations)
%    coeff:     function handle for the viscosity coefficient (mu in the equation)
%    Cpen:      a penalization term
%   
% OUTPUT:
%
%     N_mat:       the computed matrix, to be added in the left hand-side
%     N_rhs:       the computed right hand-side
%     vel:         assigned value to the normal degrees of freedom
%     normal_dofs: global numbering of the normal basis functions
%
% Copyright (C) 2014 Adriano Cortes, Rafael Vazquez
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

function [A, rhs, u, dofs] = sp_weak_drchlt_bc (space, msh, geometry, der2, bnd_sides, bnd_func, coeff, Cpen)

  A = spalloc (space.ndof, space.ndof, 3*space.ndof);
  rhs = zeros (space.ndof, 1);

  ndim = numel (msh.qn);

% Compute the matrices to impose the tangential boundary condition weakly
  for iside = bnd_sides

    if (ndim == 2)
      g_func = @(x, y) bnd_func (x, y, iside);
    elseif (ndim == 3)
      g_func = @(x, y, z) bnd_func (x, y, z, iside);
    end
    msh_side = msh_eval_boundary_side (msh, iside);

    if (ndim == 2)
      if (iside == 1)
        msh_aux = msh_2d ({msh.breaks{1}(1:2) msh.breaks{2}}, {msh.breaks{1}(1), msh.qn{2}}, {1, msh.qw{2}}, geometry, 'der2', der2);
      elseif (iside == 2)
        msh_aux = msh_2d ({msh.breaks{1}(end-1:end) msh.breaks{2}}, {msh.breaks{1}(end), msh.qn{2}}, {1, msh.qw{2}}, geometry, 'der2', der2);
      elseif (iside == 3)
        msh_aux = msh_2d ({msh.breaks{1} msh.breaks{2}(1:2)}, {msh.qn{1}, msh.breaks{2}(1)}, {msh.qw{1}, 1}, geometry, 'der2', der2);
      elseif (iside == 4)
        msh_aux = msh_2d ({msh.breaks{1} msh.breaks{2}(end-1:end)}, {msh.qn{1}, msh.breaks{2}(end)}, {msh.qw{1}, 1}, geometry, 'der2', der2);
      end
    elseif (ndim == 3)
      error ('Weak Dirichlet boundary conditions not implemented in 3D yet')
    end

    sp_bnd = space.constructor (msh_aux);
    sp_bnd = struct (sp_precompute (sp_bnd, msh_aux));
    
    for idim = 1:ndim
      x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
    end
    
    coeff_at_qnodes = coeff (x{:});

    % Since trial and test spaces are the same, we can use B'
    B = op_gradv_n_u (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);

    g_times_coeff = bsxfun (@times, g_func(x{:}), ...
         reshape(coeff_at_qnodes,[1, msh_side.nqn, msh_side.nel]));
    gradv_n_g = op_gradv_n_f (sp_bnd, msh_side, g_times_coeff);

keyboard
    coeff_at_qnodes =  coeff_at_qnodes * Cpen ./ msh_side.charlen;
    C = op_u_v (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);

    g_times_coeff = bsxfun (@times, g_func(x{:}), ...
         reshape(coeff_at_qnodes,[1, msh_side.nqn, msh_side.nel]));
    g_cdot_v = op_f_v (sp_bnd, msh_side, g_times_coeff);

    A = A + (B + B' - C);
    rhs = rhs - gradv_n_g + g_cdot_v;
  end
  
  
% The normal boundary condition is imposed strongly
  M = spalloc (space.ndof, space.ndof, 3*space.ndof);
  rhs2 = zeros (space.ndof, 1);

  normal_dofs = [];
  for iside = bnd_sides  
    msh_side = msh_eval_boundary_side (msh, iside);
    sp_side = sp_eval_boundary_side (space, msh_side);

    ind = floor ((iside+1)/2); % ind = [1 1 2 2];
    normal_dofs = union (normal_dofs, space.boundary(iside).comp_dofs{ind});

    for idim = 1:ndim
      x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
    end
    if (ndim == 2)
      g_func = @(x, y) bnd_func (x, y, iside);
    elseif (ndim == 3)
      g_func = @(x, y, z) bnd_func (x, y, z, iside);
    end
    g = g_func(x{:});
    
    M_loc = op_udotn_vdotn (sp_side, sp_side, msh_side, ones(size(x{1})));
    rhs_loc = op_fdotn_vdotn (sp_side, msh_side, g);
    
    M(sp_side.dofs, sp_side.dofs) = M(sp_side.dofs, sp_side.dofs) + M_loc;
    rhs2(sp_side.dofs) = rhs2(sp_side.dofs) + rhs_loc;
  end
  
  u = M(normal_dofs, normal_dofs) \ rhs2(normal_dofs);
  dofs = normal_dofs;
  
end
