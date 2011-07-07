% OP_CURLV_P: assemble the matrix B = [b(i,j)], b(i,j) = (coeff p_i, curl v_j).
%
%   mat = op_curlv_p (spv, spp, msh, coeff);
%   [rows, cols, values] = op_curlv_p (spv, spp, msh, coeff);
%
% INPUT:
%   
%  spv:   structure representing the space of vectorial trial functions (see sp_vector_2d_curl_transform/sp_evaluate_col)
%  spp:   structure representing the space of scalar test functions  (see sp_bspline_2d_3forms/sp_evaluate_col)
%  msh:   structure containing the domain partition and the quadrature rule (see msh_2d/msh_evaluate_col)
%  coeff: physical parameter
%
% OUTPUT:
%
%  mat:    assembled mass matrix
%  rows:   row indices of the nonzero entries
%  cols:   column indices of the nonzero entries
%  values: values of the nonzero entries
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

function varargout = op_curlv_p (spv, spp, msh, coeff)
  
  rows = zeros (msh.nel * spp.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spp.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spp.nsh_max * spv.nsh_max, 1);

  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
      jacdet_weights = msh.jacdet(:, iel) .* ...
                       msh.quad_weights(:, iel) .* coeff(:, iel);
      jacdet_weights = repmat (jacdet_weights, [1, spv.nsh(iel)]);

      shapep_iel = reshape (spp.shape_functions(:, 1:spp.nsh(iel), iel), ...
                                                  msh.nqn, spp.nsh(iel));
      curlv_iel = reshape (spv.shape_function_curls(:, 1:spv.nsh(iel), iel), ...
                           msh.nqn, spv.nsh(iel));

      for idof = 1:spp.nsh(iel)
        ishp = repmat (shapep_iel(:, idof), [1, spv.nsh(iel)]);

        rows(ncounter+(1:spv.nsh(iel))) = spp.connectivity(idof, iel);
        cols(ncounter+(1:spv.nsh(iel))) = spv.connectivity(1:spv.nsh(iel), iel);

        values(ncounter+(1:spv.nsh(iel))) = ...
               sum (jacdet_weights .* ishp .* curlv_iel, 1);

        ncounter = ncounter + spv.nsh(iel);
      end
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_curlv_p: singular map in element number %d', iel)
    end
  end

  if (nargout == 1)
    varargout{1} = sparse (rows, cols, values, spp.ndof, spv.ndof);
  elseif (nargout == 3)
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = values;
  else
    error ('op_curlv_p: wrong number of output arguments')
  end

end

