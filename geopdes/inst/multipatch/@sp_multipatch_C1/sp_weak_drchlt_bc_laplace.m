% SP_WEAK_DRCHLT_BC_LAPLACE: compute the matrix and right hand-side to impose
%  the Dirichlet boundary conditions in weak form for Laplace (Poisson) problem.
%
% The code computes the following terms in the left hand-side
% 
%  - \int_{Gamma_D} mu*{du/dn v - dv/dn u + (Cpen /  he) * (u v)}
%
% and in the right hand-side
%
%  - \int_{Gamma_D} mu*{dv/dn g + (Cpen / he) * (v g)}
%
% with u the trial function, v the test function, he the normal characteristic length, 
%  and g the boundary condition to be imposed.
%
%
%   [N_mat, N_rhs] = sp_weak_drchlt_bc_laplace  (space, msh, refs, bnd_func, coeff, Cpen)
%
% INPUTS:
%
%  space_v:    object for the multipatch space (see sp_multipatch_C1). 
%  msh:        object for the multipatch mesh (see msh_multipatch)
%  refs:       boundary sides on which the Dirichlet condition is imposed
%  bnd_func:   the condition to be imposed (g in the equations)
%  coeff:      function handle for the viscosity coefficient (mu in the equation)
%  Cpen:       a penalization term
%   
% OUTPUT:
%
%   N_mat:     the computed matrix, to be added in the left hand-side
%   N_rhs:     the computed right hand-side
%
% Copyright (C) 2014 Adriano Cortes, Rafael Vazquez
% Copyright (C) 2015, 2022 Rafael Vazquez
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

function [A, rhs] = sp_weak_drchlt_bc_laplace (space, msh, refs, bnd_func, coeff, Cpen, varargin)

  if (nargin < 6 || isempty (Cpen))
    Cpen = 10 * (min (space.sp_patch{1}.degree) + 1); 
  end

  A = spalloc (space.ndof, space.ndof, 3*space.ndof);
  rhs = zeros (space.ndof, 1);

% Compute the matrices to impose the tangential boundary condition weakly
  for iref = refs
    for bnd_side = 1:msh.boundaries(iref).nsides
      iptc = msh.boundaries(iref).patches(bnd_side);
      iside = msh.boundaries(iref).faces(bnd_side);

      msh_side = msh_eval_boundary_side (msh.msh_patch{iptc}, iside);
      msh_side_from_interior = msh_boundary_side_from_interior (msh.msh_patch{iptc}, iside);

      sp_bnd = space.sp_patch{iptc}.constructor (msh_side_from_interior);
      sp_bnd = sp_precompute (sp_bnd, msh_side_from_interior, 'value', true, 'gradient', true);

      for idim = 1:msh.rdim
        x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
      end

      coeff_at_qnodes = coeff (x{:});

      % Since trial and test spaces are the same, we can use B'
      B = op_gradv_n_u (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);

      g_times_coeff = bnd_func(x{:}, iref) .* coeff_at_qnodes;
      gradv_n_g = op_gradv_n_f (sp_bnd, msh_side, g_times_coeff);

      coeff_at_qnodes =  coeff_at_qnodes * Cpen ./ msh_side.charlen;
      C = op_u_v (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);

      g_times_coeff = bnd_func(x{:}, iref) .* coeff_at_qnodes;
      g_cdot_v = op_f_v (sp_bnd, msh_side, g_times_coeff);

      [Cpatch, Cpatch_cols] = sp_compute_Cpatch (space, iptc);
      A(Cpatch_cols,Cpatch_cols) = A(Cpatch_cols,Cpatch_cols) + ...
        Cpatch.' * (B + B.' - C) * Cpatch;
      rhs(Cpatch_cols) = rhs(Cpatch_cols) + Cpatch.' * (-gradv_n_g + g_cdot_v);
%       dofs = space.gnum{iptc};
%       A(dofs,dofs) = A(dofs,dofs) + (B + B' - C);
%       rhs(dofs) = rhs(dofs) - gradv_n_g + g_cdot_v;
    end
  end

end
