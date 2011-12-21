% OP_CURLU_CURLV_3D: assemble the matrix A = [a(i,j)], a(i,j) = (coeff curl u_j, curl v_i), with vector-valued curl.
%
%   mat = op_curlu_curlv_3d (spu, spv, msh, epsilon);
%   [rows, cols, values] = op_curlu_curlv_3d (spu, spv, msh, epsilon);
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_vector_3d_curl_transform/sp_evaluate_col)
%   spv:   structure representing the space of test functions  (see sp_vector_3d_curl_transform/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_3d/msh_evaluate_col)
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


function varargout = op_curlu_curlv_3d (spu, spv, msh, coeff)
  
  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
      jacdet_weights = reshape (msh.jacdet(:, iel) .* ...
                       msh.quad_weights(:, iel) .* coeff(:, iel), 1, msh.nqn);

      curlv_iel = reshape (spv.shape_function_curls(:, :, 1:spv.nsh(iel), iel), ...
                            3, msh.nqn, spv.nsh(iel));
      curlu_iel = reshape (spu.shape_function_curls(:, :, 1:spu.nsh(iel), iel), ...
                            3, msh.nqn, spu.nsh(iel));

      curlv_times_jw = bsxfun (@times, jacdet_weights, curlv_iel);
      for idof = 1:spv.nsh(iel)
        rows(ncounter+(1:spu.nsh(iel))) = spv.connectivity(idof, iel);
        cols(ncounter+(1:spu.nsh(iel))) = spu.connectivity(1:spu.nsh(iel), iel);

        aux_val = bsxfun (@times, curlv_times_jw(:,:,idof), curlu_iel);
        values(ncounter+(1:spu.nsh(iel))) = sum (sum (aux_val, 2), 1);
        ncounter = ncounter + spu.nsh(iel);
      end
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_curlu_curlv_3d: singular map in element number %d', iel)
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
    error ('op_curlu_curlv_3d: wrong number of output arguments')
  end

end

%% COPY OF THE FIRST VERSION OF THE FUNCTION (MORE UNDERSTANDABLE)
% 
% function mat = op_curlu_curlv_3d (spu, spv, msh, coeff)
%   
%   mat = spalloc(spv.ndof, spu.ndof, 1);
%   for iel = 1:msh.nel
%     if (all (msh.jacdet(:,iel)))
%       mat_loc = zeros (spv.nsh(iel), spu.nsh(iel));
%       for idof = 1:spv.nsh(iel)
%         ishc = reshape(spv.shape_function_curls(:,:,idof,iel), 3, []);
%         for jdof = 1:spu.nsh(iel)
%           jshc = reshape (spu.shape_function_curls(:,:,jdof,iel), 3, []);
% % The cycle on the quadrature points is vectorized
%           %for inode = 1:msh.nqn
%             mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
%               sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .*...
%               sum (ishc .* jshc, 1).' .* coeff (:, iel));
%           %end  
%         end
%       end
%       mat(spv.connectivity(:,iel), spu.connectivity(:,iel)) = ...
%         mat(spv.connectivity(:,iel), spu.connectivity(:,iel)) + mat_loc;
%     else
%       warning ('geopdes:jacdet_zero_at_quad_node', 'op_curlu_curlv_3d: singular map in element number %d', iel)
%     end
%   end
% 
% end