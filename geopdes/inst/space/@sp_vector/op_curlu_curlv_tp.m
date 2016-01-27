% OP_CURLU_CURLV_TP: assemble the matrix A = [a(i,j)], a(i,j) = (epsilon curl u_j, curl v_i), exploiting the tensor product structure.
%
%   mat = op_curlu_curlv_tp (spu, spv, msh, [epsilon]);
%   [rows, cols, values] = op_curlu_curlv_tp (spu, spv, msh, [epsilon]);
%
% The same function works for 2d (scalar-valued curl) and 3d problems (vector-valued curl).
%
% INPUT:
%
%   spu:     object that defines the space of trial functions (see sp_vector)
%   spv:     object that defines the space of test functions (see sp_vector)
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
% Copyright (C) 2011, 2015, 2016 Rafael Vazquez
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

function varargout = op_curlu_curlv_tp (space1, space2, msh, coeff)

  A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);

  if (msh.ndim == 2)
    for iel = 1:msh.nel_dir(1)
      msh_col = msh_evaluate_col (msh, iel);
      sp1_col = sp_evaluate_col (space1, msh_col, 'value', false, 'curl', true);
      sp2_col = sp_evaluate_col (space2, msh_col, 'value', false, 'curl', true);

      if (nargin == 4)
        for idim = 1:msh.rdim
          x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
        end
        coeffs = coeff (x{:});
      else
        coeffs = ones (msh_col.nqn, msh_col.nel);
      end

      A = A + op_curlu_curlv_2d (sp1_col, sp2_col, msh_col, coeffs);
    end
  elseif (msh.ndim == 3)
    for iel = 1:msh.nel_dir(1)
      msh_col = msh_evaluate_col (msh, iel);
      sp1_col = sp_evaluate_col (space1, msh_col, 'value', false, 'curl', true);
      sp2_col = sp_evaluate_col (space2, msh_col, 'value', false, 'curl', true);

      if (nargin == 4)
        for idim = 1:msh.rdim
          x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
        end
        coeffs = coeff (x{:});
      else
        coeffs = ones (msh_col.nqn, msh_col.nel);
      end

      A = A + op_curlu_curlv_3d (sp1_col, sp2_col, msh_col, coeffs);
    end
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
