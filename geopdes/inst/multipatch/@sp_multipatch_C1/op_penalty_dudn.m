% OP_PENALTY_DUDN: matrix and right-hand side to compute penalty terms to impose du/dn = f.
%  It computes the terms of the form Cp*(du/dn, dv/dn) and Cp*(f, dv/dn), 
%  using the same space for trial and test functions.
%
%   [mat, rhs] = op_penalty_dudn (space, msh, sides, Cpen, [coeff]);
%
% INPUT:
%
%  space: object representing the space of trial functions (see sp_multipatch_C1)
%  msh:   object defining the domain partition and the quadrature rule (see msh_multipatch)
%  sides: boundary sides on which to compute the integrals
%  Cpen:  penalization parameter, Cp in the equation above
%  coeff: function handle to compute the Neumann condition. If empty, the returned rhs will be zero.
%
% OUTPUT:
%
%  mat:    assembled mass matrix
%  rhs:    assembled right-hand side
% 
% Copyright (C) 2023, 2024 Michele Torre, Rafael Vazquez
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

function [P, rhs] = op_penalty_dudn (space, msh, sides, Cpen, coeff)

  P =  spalloc (space.ndof, space.ndof, space.ndof);
  rhs = zeros(space.ndof, 1);

  for iref = sides
    for bnd_side = 1:msh.boundaries(iref).nsides
      iptc = msh.boundaries(iref).patches(bnd_side);
      iside = msh.boundaries(iref).faces(bnd_side);

      msh_side = msh_eval_boundary_side (msh.msh_patch{iptc}, iside);
      msh_side_from_interior = msh_boundary_side_from_interior (msh.msh_patch{iptc}, iside);

      sp_bnd = space.sp_patch{iptc}.constructor (msh_side_from_interior);
      sp_bnd = sp_precompute (sp_bnd, msh_side_from_interior, 'value', false, 'gradient', true);

      coe_side = Cpen .* msh_side.charlen; 
      tmp = op_gradu_n_gradv_n(sp_bnd, sp_bnd, msh_side, coe_side);
      
      [Cpatch, Cpatch_cols] = sp_compute_Cpatch (space, iptc);
      P(Cpatch_cols,Cpatch_cols) = P(Cpatch_cols,Cpatch_cols) + ...
        Cpatch.' * (tmp) * Cpatch;

      if (nargin == 5)
        x = cell (msh.rdim, 1);
        for idim = 1:msh.rdim
          x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
        end
        rhs(Cpatch_cols) = rhs(Cpatch_cols) + ...
          Cpatch.' *op_gradv_n_f (sp_bnd, msh_side, Cpen*coeff(x{:}));
      end
      
    end
  end

end
