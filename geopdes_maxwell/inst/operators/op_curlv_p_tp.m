% OP_CURLV_P: assemble the matrix B = [b(i,j)], b(i,j) = (coeff p_i, curl v_j), exploiting the tensor product structure.
%
%   mat = op_curlv_p_tp (spv, spp, msh, coeff);
%   [rows, cols, values] = op_curlv_p (spv, spp, msh, coeff);

% INPUT:
%
%   spv:     object that defines the space of trial functions (see sp_vector_2d_curl_transform)
%   spp:     object that defines the space of the multiplier (see sp_bspline_2d_3forms)
%   msh:     object that defines the domain partition and the quadrature rule (see msh_2d)
%   epsilon: function handle to compute some physical coefficient
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
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

function varargout = op_curlv_p_tp (spv, spp, msh, coeff)

  A = spalloc (spp.ndof, spv.ndof, 3*spv.ndof);

  ndim = numel (msh.qn);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    spv_col = sp_evaluate_col (spv, msh_col, 'value', false, 'curl', true);
    spp_col = sp_evaluate_col (spp, msh_col);

    for idim = 1:ndim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end

    A = A + op_curlv_p (spv_col, spp_col, msh_col, coeff (x{:}));
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
