% OP_DIV_V_Q_TP: assemble the matrix B = [b(i,j)], b(i,j) = (q_i, div v_j), exploiting the tensor product structure.
%
%   mat = op_div_v_q_tp (spv, spq, msh);
%   [rows, cols, values] = op_div_v_q_tp (spv, spq, msh);
%
% INPUT: 
%
%   spv:     object representing the vector-valued space of trial functions for the velocity (see sp_vector)
%   spq:     object representing the scalar-valued space of test functions for the pressure (see sp_scalar)
%   msh:     object that defines the domain partition and the quadrature rule (see msh_cartesian)
%
% OUTPUT: 
%
%   mat:    assembled matrix 
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011 Andrea Bressan
% Copyright (C) 2011, 2017 Rafael Vazquez
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

  for idim = 1:msh.ndim
    size2 = size (spq.sp_univ(idim).connectivity);
    for icomp = 1:spv.ncomp_param
      size1 = size (spv.scalar_spaces{icomp}.sp_univ(idim).connectivity);
      if (size1(2) ~= size2(2) || size1(2) ~= msh.nel_dir(idim))
        error ('One of the discrete spaces is not associated to the mesh')
      end
    end
  end

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
