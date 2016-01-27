% OP_V_GRADP_TP: assemble the matrix B = [b(i,j)], b(i,j) = (epsilon grad p_i, v_j), exploiting the tensor product structure.
%
%   mat = op_v_gradp_tp (spv, spp, msh, [epsilon]);
%   [rows, cols, values] = op_v_gradp_tp (spv, spp, msh, [epsilon]);
%
% INPUT:
%
%   spv:     object that defines the vector-valued space of trial functions (see sp_vector)
%   spp:     object that defines the scalar-valued space of the multiplier (see sp_scalar)
%   msh:     object that defines the domain partition and the quadrature rule (see msh_cartesian)
%   epsilon: function handle to compute some physical coefficient (optional)
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2011, 2016 Rafael Vazquez
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

function varargout = op_v_gradp_tp (spv, spp, msh, coeff)

  A = spalloc (spp.ndof, spv.ndof, 3*spv.ndof);

  ndim = numel (msh.qn);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    spv_col = sp_evaluate_col (spv, msh_col);
    spp_col = sp_evaluate_col (spp, msh_col, 'value', false, 'gradient', true);

    if (nargin == 4)
      for idim = 1:msh.rdim
        x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
      end
      coeffs = coeff (x{:});
    else
      coeffs = ones (msh_col.nqn, msh_col.nel);
    end

    A = A + op_v_gradp (spv_col, spp_col, msh_col, coeffs);
  end

  if (nargout == 1)
    varargout{1} = A;
  elseif (nargout == 3)
    [rows, cols, vals] = find (A);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end

end
