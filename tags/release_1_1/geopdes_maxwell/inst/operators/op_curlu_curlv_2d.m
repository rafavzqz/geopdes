% OP_CURLU_CURLV_2D: assemble the matrix A = [a(i,j)], a(i,j) = (coeff curl u_j, curl v_i), with scalar-valued curl.
%
%   mat = op_curlu_curlv_2d (spu, spv, msh, epsilon);
%   [rows, cols, values] = op_curlu_curlv_2d (spu, spv, msh, epsilon);
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_bsp_hcurl_2d_phys)
%   spv:   structure representing the space of test functions  (see sp_bsp_hcurl_2d_phys)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%   coeff: physical parameter
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


function varargout = op_curlu_curlv_2d (spu, spv, msh, coeff)
  
  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
      jacdet_weights = msh.jacdet(:,iel) .* ...
                       msh.quad_weights(:, iel) .* coeff (:, iel);
      for idof = 1:spv.nsh(iel)
        ishc = squeeze (spv.shape_function_curls(:, idof, iel));
        for jdof = 1:spu.nsh(iel)
          ncounter = ncounter + 1;
          rows(ncounter) = spv.connectivity(idof, iel);
          cols(ncounter) = spu.connectivity(jdof, iel);

          jshc = squeeze (spu.shape_function_curls(:, jdof, iel));

          values(ncounter) = sum (jacdet_weights .* ishc .* jshc);
        end
      end
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_curlu_curlv: singular map in element number %d', iel)
    end
  end

  if (nargout == 1)
    varargout{1} = sparse (rows, cols, values, spv.ndof, spu.ndof);
  elseif (nargout == 3)
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = values;
  else
    error ('op_curlu_curlv_2d: wrong number of output arguments')
  end

end

