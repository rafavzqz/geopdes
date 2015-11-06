% OP_DIV_V_Q_TP: assemble the matrix B = [b(i,j)], b(i,j) = (q_i, div v_j), exploiting the tensor product structure.
%
%   mat = op_div_v_q_tp (spv, spq, msh);
%   [rows, cols, values] = op_div_v_q_tp (spv, spq, msh);
%
% INPUT: 
%
%   spv:     object representing the space of trial functions for the velocity (see sp_vector_2d, sp_vector_2d_piola)
%   spq:     object representing the space of test functions for the pressure (see sp_bspline_2d)
%   msh:     object that defines the domain partition and the quadrature rule (see msh_2d)
%
% OUTPUT: 
%
%   mat:    assembled matrix 
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011 Rafael Vazquez, Andrea Bressan
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

function varargout = op_div_v_q_tp (spv, spq, msh)

  A = spalloc (spq.ndof, spv.ndof, 3*spv.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    spv_col = sp_evaluate_col (spv, msh_col, 'divergence', true, ...
			       'value', false);
    spq_col = sp_evaluate_col (spq, msh_col);

    A = A + op_div_v_q (spv_col, spq_col, msh_col);
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
