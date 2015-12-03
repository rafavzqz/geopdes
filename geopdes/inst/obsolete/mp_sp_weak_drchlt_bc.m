% MP_SP_WEAK_DRCHLT_BC: compute the matrix and right hand-side to impose
%  the Dirichlet boundary conditions in weak form for the tangential
%  component. To be used with the 'RT' and 'NDL' spaces.
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
%   [N_mat, N_rhs] = mp_sp_weak_drchlt_bc  (space, msh, gnum, ornt, boundaries, refs, bnd_func, coeff, Cpen)
%
% INPUTS:
%
%  space:      space object (see sp_vector_div_transform)
%  msh:        mesh object (see msh_cartesian)
%  gnum:       global numbering of the degrees of freedom (see mp_interface_hdiv)
%  ornt:       global orientation of the degrees of freedom (see mp_interface_hdiv)
%  boundaries: array of structures containing the information for the boundaries (see mp_geo_load)
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
% Copyright (C) 2015 Rafael Vazquez
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

function [A, rhs] = mp_sp_weak_drchlt_bc (space, msh, gnum, ornt, boundaries, refs, bnd_func, coeff, Cpen)

  if (isa (space, 'sp_multipatch'))
    warning ('For spaces of the class SP_MULTIPATCH, using the function SP_WEAK_DRCHLT_BC inside the class')
    [A, rhs] = mp_sp_weak_drchlt_bc (space, msh, refs, bnd_func, coeff, Cpen);
    return
  end
  
  ndof = max ([gnum{:}]);

  A = spalloc (ndof, ndof, 3*ndof);
  rhs = zeros (ndof, 1);

% Compute the matrices to impose the tangential boundary condition weakly
  for iref = refs
    for bnd_side = 1:boundaries(iref).nsides
      iptc = boundaries(iref).patches(bnd_side);
      iside = boundaries(iref).faces(bnd_side);

      msh_side = msh_eval_boundary_side (msh{iptc}, iside);
      msh_side_from_interior = msh_boundary_side_from_interior (msh{iptc}, iside);

      sp_bnd = space{iptc}.constructor (msh_side_from_interior);
      sp_bnd = struct (sp_precompute (sp_bnd, msh_side_from_interior, 'value', true, 'gradient', true));

      for idim = 1:msh{iptc}.rdim
        x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
      end

      coeff_at_qnodes = coeff (x{:});

      % Since trial and test spaces are the same, we can use B'
      B = op_gradv_n_u (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);

      g_times_coeff = bsxfun (@times, bnd_func(x{:}, iref), ...
           reshape(coeff_at_qnodes,[1, msh_side.nqn, msh_side.nel]));
      gradv_n_g = op_gradv_n_f (sp_bnd, msh_side, g_times_coeff);

      coeff_at_qnodes =  coeff_at_qnodes * Cpen ./ msh_side.charlen;
      C = op_u_v (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);

      g_times_coeff = bsxfun (@times, bnd_func(x{:}, iref), ...
           reshape(coeff_at_qnodes,[1, msh_side.nqn, msh_side.nel]));
      g_cdot_v = op_f_v (sp_bnd, msh_side, g_times_coeff);

      dofs = gnum{iptc};
      ornt_matrix = spdiags (ornt{iptc}', 0, sp_bnd.ndof, sp_bnd.ndof);
      A(dofs,dofs) = A(dofs,dofs) + ornt_matrix * (B + B' - C) * ornt_matrix;
      rhs(dofs) = rhs(dofs) + (-gradv_n_g + g_cdot_v) .* ornt{iptc}';
    end
  end

end
