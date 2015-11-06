% OP_UXN_VXN_2D: assemble the matrix M = [m(i,j)], m(i,j) = (mu u_j x n, v_i x n), with n the exterior normal vector.
%
%   mat = op_uxn_vxn_2d (spu, spv, msh, coeff);
%   [rows, cols, values] = op_uxn_vxn_2d (spu, spv, msh, coeff);
%
% INPUT:
%   
%  spu:   structure representing the space of trial functions (see sp_vector_2d_curl_transform/sp_eval_boundary_side)
%  spv:   structure representing the space of test functions  (see sp_vector_2d_curl_transform/sp_eval_boundary_side)
%  msh:   structure containing the domain partition and the quadrature rule (see msh_2d/msh_eval_boundary_side)
%  coeff: physical parameter
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2009, 2010 Carlo de Falco, Rafael Vazquez
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

function varargout = op_uxn_vxn_2d (spu, spv, msh, coeff)
  
  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      jacdet_weights = msh.jacdet(:, iel) .* ...
              msh.quad_weights(:, iel) .* coeff(:, iel);

      for idof = 1:spv.nsh(iel)
        ishp = squeeze (spv.shape_functions(:, :, idof, iel));
        ishp_x_n = (ishp(1, :) .* msh.normal(2, :, iel) - ...
                    ishp(2, :) .* msh.normal(1, :, iel))';

        rows(ncounter+(1:spu.nsh(iel))) = spv.connectivity(idof, iel);
        cols(ncounter+(1:spu.nsh(iel))) = spu.connectivity(1:spu.nsh(iel), iel);

        for jdof = 1:spu.nsh(iel)
          jshp = squeeze (spu.shape_functions(:, :, jdof, iel));
          jshp_x_n = (jshp(1, :) .* msh.normal(2, :, iel) - ...
                      jshp(2, :) .* msh.normal(1, :, iel))';

          values(ncounter+jdof) = sum (jacdet_weights .* ishp_x_n .* jshp_x_n);
        end
        ncounter = ncounter + spu.nsh(iel);
      end
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_uxn_vxn_2d: singular map in element number %d', iel)
    end
  end

  if (nargout == 1)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), spv.ndof, spu.ndof);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_uxn_vxn_2d: wrong number of output arguments')
  end

end

