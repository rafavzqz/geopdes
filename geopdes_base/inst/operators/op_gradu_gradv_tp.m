% OP_GRADU_GRADV_TP: assemble the stiffness matrix A = [a(i,j)], a(i,j) = (epsilon grad u_j, grad v_i), using a tensor product space.
%
%   mat = op_gradu_gradv_tp (spu, spv, msh, epsilon);
%   [rows, cols, values] = op_gradu_gradv_tp (spu, spv, msh, epsilon);
%
% INPUT:
%
%   spu:     class representing the space of trial functions (see sp_bspline_2d)
%   spv:     class representing the space of test functions (see sp_bspline_2d)
%   msh:     structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%   epsilon: diffusion coefficient
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011, Rafael Vazquez
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

function varargout = op_gradu_gradv_tp (space1, space2, msh, coeff)

  A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);

  for iel = 1:msh.nelu
    [sp1, element_list] = sp_evaluate_col (space1, msh, iel, 'value', false);
    sp2 = sp_evaluate_col (space2, msh, iel, 'value', false);

    A = A + op_gradu_gradv (sp1, sp2, msh, coeff, element_list);
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
