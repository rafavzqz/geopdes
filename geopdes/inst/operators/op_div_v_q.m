% OP_DIV_V_Q: assemble the matrix B = [b(i,j)], b(i,j) = (q_i, div v_j).
%
%   mat = op_div_v_q (spv, spq, msh);
%   [rows, cols, values] = op_div_v_q (spv, spq, msh);
%
% INPUT: 
%
%   spv:     structure representing the space of trial functions for the velocity (see sp_vector/sp_evaluate_col)
%   spq:     structure representing the space of test functions for the pressure (see sp_scalar/sp_evaluate_col) 
%   msh:     structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%
% OUTPUT: 
%
%   mat:    assembled matrix 
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2009, 2010 Carlo de Falco
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

function varargout = op_div_v_q (spv, spq, msh)
  
  rows = zeros (msh.nel * spq.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spq.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spq.nsh_max * spv.nsh_max, 1);

  jacdet_weights = msh.jacdet .* msh.quad_weights;
  
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      divv_iel = reshape (spv.shape_function_divs(:, :, iel), msh.nqn, 1, spv.nsh_max);
      shpq_iel = reshape (spq.shape_functions(:, :, iel), msh.nqn, spq.nsh_max, 1);

      jacdet_iel = reshape (jacdet_weights(:,iel), [msh.nqn,1,1]);

      jacdet_divv = bsxfun (@times, jacdet_iel, divv_iel);
      tmp1 = bsxfun (@times, jacdet_divv, shpq_iel);
      elementary_values = reshape (sum (tmp1, 1), spq.nsh_max, spv.nsh_max);

      [rows_loc, cols_loc] = ndgrid (spq.connectivity(:,iel), spv.connectivity(:,iel));
      indices = rows_loc & cols_loc;
      rows(ncounter+(1:spv.nsh(iel)*spq.nsh(iel))) = rows_loc(indices);
      cols(ncounter+(1:spv.nsh(iel)*spq.nsh(iel))) = cols_loc(indices);
      values(ncounter+(1:spv.nsh(iel)*spq.nsh(iel))) = elementary_values(indices);
      ncounter = ncounter + spv.nsh(iel)*spq.nsh(iel);
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_div_v_q: singular map in element number %d', iel)
    end
  end

  if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), spq.ndof, spv.ndof);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_div_v_q: wrong number of output arguments')
  end

end


%% COPY OF THE FIRST VERSION OF THE FUNCTION (MORE UNDERSTANDABLE)
% 
% function mat = op_div_v_q (spv, spq, msh)
%   
%   mat = spalloc(spq.ndof, spv.ndof, 1);
%   
%   for iel = 1:msh.nel
%     if (all (msh.jacdet(:, iel)))
%       mat_loc = zeros (spq.nsh(iel), spv.nsh(iel));
%       for idof = 1:spq.nsh(iel)
%         ishp = spq.shape_functions(:, idof, iel);
%         for jdof = 1:spv.nsh(iel) 
%           jshd = spv.shape_function_divs(:, jdof, iel);
% % The cycle on the quadrature points is vectorized
%           %for inode = 1:msh.nqn
%           mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
%             sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* ishp .* jshd);
%           %end  
%         end
%       end
%       mat(spq.connectivity(:, iel), spv.connectivity(:, iel)) = ...
%         mat(spq.connectivity(:, iel), spv.connectivity(:, iel)) + mat_loc;
%     else
%       warning ('geopdes:jacdet_zero_at_quad_node', 'op_div_v_q: singular map in element number %d', iel)
%     end
%   end
% 
% end
