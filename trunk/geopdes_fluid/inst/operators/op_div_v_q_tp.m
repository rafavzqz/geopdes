% OP_DIV_V_Q_TP: assemble the matrix B = [b(i,j)], b(i,j) = (q_i, div v_j).
%
%   mat = op_div_v_q (spv, spq, msh);
%   [rows, cols, values] = op_div_v_q (spv, spq, msh);
%
% INPUT: 
%
%   spv:     structure representing the space of trial functions for the velocity (see sp_bspline_2d_phys) 
%   spq:     structure representing the space of test functions for the pressure (see sp_bspline_2d_phys) 
%   msh:     structure containing the domain partition and the quadrature rule (see msh_push_forward_2d) 
%
% OUTPUT: 
%
%   mat:    assembled matrix 
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011 Rafael Vazquez
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

function varargout = op_div_v_q_tp (space1, space2, msh)

  A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);

  for iel = 1:msh.nelu
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space1, msh_col, 'divergence', true);
    sp2_col = sp_evaluate_col (space2, msh_col, 'gradient', false);

    A = A + op_div_v_q (sp1_col, sp2_col, msh_col);
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
