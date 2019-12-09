% OP_KL_SHELLS_TP: assemble the Kirchhoff-Love shell stiffness matrix, exploiting the tensor product structure.
%
%   mat = op_KL_shells_tp (spu, spv, msh, E_coeff, nu_coeff, t_coeff);
%   [rows, cols, values] = op_KL_shells_tp (spu, spv, msh, E_coeff, nu_coeff, t_coeff);
%
% INPUT:
%
%  spu:   object representing the space of trial functions (see sp_vector)
%  spv:   object representing the space of test functions (see sp_vector)
%  msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%  E_coeff:  function handle to compute the Young's modulus
%  nu_coeff: function handle to compute the Poisson's ratio
%  t_coeff:  thickness of the shell, scalar value

% OUTPUT:
%
%  mat:    assembled stiffness matrix
%  rows:   row indices of the nonzero entries
%  cols:   column indices of the nonzero entries
%  values: values of the nonzero entries
% 
% Copyright (C) 2018, 2019 Pablo Antolin, Luca Coradello
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

function varargout = op_KL_shells_tp (space_u, space_v, msh, E_coeff, nu_coeff, t_coeff)

  A = spalloc (space_v.ndof, space_u.ndof, 3*space_u.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col_u = sp_evaluate_col_param (space_u, msh_col, 'value', true, 'gradient', true, 'hessian', true);
    sp_col_v = sp_evaluate_col_param (space_v, msh_col, 'value', true, 'gradient', true, 'hessian', true);
    
    if (nargin == 6)
      for idim = 1:msh.rdim
        x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
      end
      E_coeffs = E_coeff (x{:});
      nu_coeffs = nu_coeff (x{:});
    else
      error ('geopdes:op_KL_shells_tp: invalid input')
    end

    A = A + op_KL_shells (sp_col_u, sp_col_v, msh_col, E_coeffs, nu_coeffs, t_coeff);
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
