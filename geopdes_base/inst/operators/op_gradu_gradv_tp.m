% OP_GRADU_GRADV_TP: assemble the stiffness matrix A = [a(i,j)], a(i,j) = (epsilon grad u_j, grad v_i), exploiting the tensor product structure.
%
%   mat = op_gradu_gradv_tp (spu, spv, msh, epsilon);
%   [rows, cols, values] = op_gradu_gradv_tp (spu, spv, msh, epsilon);
%   mat = op_gradu_gradv_tp (spu, spv, msh, eps_sp, eps_dofs);
%   [rows, cols, values] = op_gradu_gradv_tp (spu, spv, msh, eps_sp, eps_dofs);
%
% INPUT:
%
%   spu:      class representing the space of trial functions (see sp_bspline_2d)
%   spv:      class representing the space of test functions (see sp_bspline_2d)
%   msh:      class defining the domain partition and the quadrature rule (see msh_2d)
%   epsilon:  function handle to compute the diffusion coefficient
%   eps_sp:   function space of the source function
%   eps_dofs: weitghts of the dofs in the source function
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2011 Rafael Vazquez
% Copyright (C) 2011, 2013 Carlo de Falco
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

function varargout = op_gradu_gradv_tp (space1, space2, msh, varargin)

  if (numel (varargin) < 2)
    fun_coef = true;
    coeff = varargin{1};
  elseif (numel (varargin) == 2)
    fun_coef = false;      
    [coeff_sp, coeff_dofs] = deal (varargin{:});
  else
    error ('op_gradu_gradv_tp: wrong number of input parameters');
  end


  A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);

  ndim = numel (msh.qn);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space1, msh_col, 'value', false, 'gradient', true);
    sp2_col = sp_evaluate_col (space2, msh_col, 'value', false, 'gradient', true);

    if (fun_coef)
      for idim = 1:ndim
        x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
      end

      A = A + op_gradu_gradv (sp1_col, sp2_col, msh_col, coeff (x{:}));

    else
      
      A = A + op_gradu_gradv (sp1_col, sp2_col, msh_col, ...
                              sp_eval (eps_dofs, eps_sp, msh_col));
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
