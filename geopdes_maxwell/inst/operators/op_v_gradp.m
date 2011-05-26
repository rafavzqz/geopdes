% OP_V_GRADP: assemble the matrix B = [b(i,j)], b(i,j) = (epsilon grad p_i, v_j).
%
%   mat = op_v_gradp (spv, spp, msh, epsilon);
%   [rows, cols, values] = op_v_gradp (spv, spp, msh, epsilon);
%
% INPUT:
%    
%   spv:     structure representing the space of vectorial trial functions  (see sp_bsp_hcurl_2d_phys)
%   spp:     structure representing the space of scalar test functions (see sp_bspline_2d_phys)
%   msh:     structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%   epsilon: physical parameter
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

function varargout = op_v_gradp (spv, spp, msh, coeff)

  ndir = size (spp.shape_function_gradients, 1);

  rows = zeros (msh.nel * spp.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spp.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spp.nsh_max * spv.nsh_max, 1);

  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      jacdet_weights = msh.jacdet(:, iel) .* ...
                       msh.quad_weights(:, iel) .* coeff(:, iel);
      jacdet_weights = repmat (jacdet_weights, [1, spv.nsh(iel)]);

      gradp_iel = reshape (spp.shape_function_gradients(:, :, 1:spp.nsh(iel), iel), ...
                                     ndir, msh.nqn, spp.nsh(iel));
      shpv_iel = reshape (spv.shape_functions(:, :, 1:spv.nsh(iel), iel), ...
                                     ndir, msh.nqn, spv.nsh(iel));

      for idof = 1:spp.nsh(iel)
        ishg  = repmat (gradp_iel(:, :, idof), [1, 1, spv.nsh(iel)]);

        rows(ncounter+(1:spv.nsh(iel))) = spp.connectivity(idof, iel);
        cols(ncounter+(1:spv.nsh(iel))) = spv.connectivity(1:spv.nsh(iel), iel);

        values(ncounter+(1:spv.nsh(iel))) = sum (jacdet_weights .* ...
          reshape (sum (ishg .* shpv_iel, 1), msh.nqn, spv.nsh(iel)), 1);
        ncounter = ncounter + spv.nsh(iel);
      end
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_v_gradp: singular map in element number %d', iel)
    end
  end

  if (nargout == 1)
    varargout{1} = sparse (rows, cols, values, spp.ndof, spv.ndof);
  elseif (nargout == 3)
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = values;
  else
    error ('op_v_gradp: wrong number of output arguments')
  end

end

