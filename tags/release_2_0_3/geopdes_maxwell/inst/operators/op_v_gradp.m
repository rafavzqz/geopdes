% OP_V_GRADP: assemble the matrix B = [b(i,j)], b(i,j) = (epsilon grad p_i, v_j).
%
%   mat = op_v_gradp (spv, spp, msh, epsilon);
%   [rows, cols, values] = op_v_gradp (spv, spp, msh, epsilon);
%
% INPUT:
%    
%   spv:     structure representing the space of vectorial trial functions  (see sp_vector_2d_curl_transform/sp_evaluate_col)
%   spp:     structure representing the space of scalar test functions (see sp_bspline_2d/sp_evaluate_col)
%   msh:     structure containing the domain partition and the quadrature rule (see msh_2d/msh_evaluate_col)
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
      jacdet_weights = reshape (msh.jacdet(:,iel) .* ...
                         msh.quad_weights(:, iel) .* coeff(:,iel), 1, msh.nqn);

      gradp_iel = reshape (spp.shape_function_gradients(:, :, 1:spp.nsh(iel), iel), ...
                                     ndir, msh.nqn, spp.nsh(iel));
      shpv_iel = reshape (spv.shape_functions(:, :, 1:spv.nsh(iel), iel), ...
                                     ndir, msh.nqn, spv.nsh(iel));

      gradp_times_jw = bsxfun (@times, jacdet_weights, gradp_iel);
      for idof = 1:spp.nsh(iel)
        rows(ncounter+(1:spv.nsh(iel))) = spp.connectivity(idof, iel);
        cols(ncounter+(1:spv.nsh(iel))) = spv.connectivity(1:spv.nsh(iel), iel);

        aux_val = bsxfun (@times, gradp_times_jw(:,:,idof), shpv_iel);
        values(ncounter+(1:spv.nsh(iel))) = sum (sum (aux_val, 2), 1);
        ncounter = ncounter + spv.nsh(iel);
      end
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_v_gradp: singular map in element number %d', iel)
    end
  end

  if (nargout == 1)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), spp.ndof, spv.ndof);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_v_gradp: wrong number of output arguments')
  end

end


%% COPY OF THE FIRST VERSION OF THE FUNCTION (MORE UNDERSTANDABLE)
% function mat = op_v_gradp (spv, spp, msh, coeff)
% 
%   mat = spalloc(spp.ndof, spv.ndof, 1);
%   ndir = size (spp.shape_function_gradients, 1);
%   for iel = 1:msh.nel
%     if (all (msh.jacdet(:,iel)))
%       mat_loc = zeros (spp.nsh(iel), spv.nsh(iel));
%       for idof = 1:spp.nsh(iel)
%         ishg = reshape(spp.shape_function_gradients(:,:,idof,iel),ndir,[]);
%         for jdof = 1:spv.nsh(iel)
%           jshg = reshape(spv.shape_functions(:,:,jdof,iel),ndir,[]);
% % The cycle on the quadrature points is vectorized
%           %for inode = 1:msh.nqn
%             mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
%              sum (msh.jacdet(:, iel) .* msh.quad_weights(:, iel) .* ...
%              sum (ishg .* jshg, 1).' .* coeff(:, iel));
%           %end  
%         end
%       end
%       mat(spp.connectivity(:, iel), spv.connectivity(:, iel)) = ...
%         mat(spp.connectivity(:, iel), spv.connectivity(:, iel)) + mat_loc;
%     else
%       warning ('geopdes:jacdet_zero_at_quad_node', 'op_v_gradp: singular map in element number %d', iel)
%     end
%   end
% 
% end