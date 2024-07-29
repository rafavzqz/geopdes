% OP_PENALTY_DUDN: matrix and right-hand side to compute penalty terms to impose du/dn = f.
%  It computes the terms of the form Cp*(du/dn, dv/dn) and Cp*(f, dv/dn), 
%  using the same space for trial and test functions.
%
%   [mat, rhs] = op_penalty_dudn (space, msh, sides, Cpen, [coeff]);
%
% INPUT:
%
%  space: object representing the space of trial functions (see sp_scalar)
%  msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
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

  for iside = 1:numel(sides)
    msh_side = msh_eval_boundary_side (msh, iside);
    msh_side_int = msh_boundary_side_from_interior (msh, iside);
    sp_side = space.constructor (msh_side_int);
    sp_side = sp_precompute (sp_side , msh_side_int , 'gradient', true);

    coe_side = Cpen .* msh_side.charlen; 
    mass_pen = op_gradu_n_gradv_n (sp_side, sp_side, msh_side, coe_side);

    if (nargin == 5)
      x = cell (msh.rdim, 1);
      for idim = 1:msh.rdim
        x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
      end
      rhs = rhs + op_gradv_n_f (sp_side, msh_side, Cpen*coeff(x{:}));
    end

    P = P + mass_pen;
  end

end
