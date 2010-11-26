% OP_DIV_V_Q: assemble the matrix B = [b(i,j)], b(i,j) = (q_i, div v_j).
%
%   mat = op_div_v_q (spv, spq, msh);
%
% INPUT: 
%
%   spv:     structure representing the space of trial functions for the velocity (see sp_bspline_2d_phys) 
%   spq:     structure representing the space of test functions for the pressure (see sp_bspline_2d_phys) 
%   msh:     structure containing the domain partition and the quadrature rule (see msh_push_forward_2d) 
%
% OUTPUT: 
%
%   mat: assembled matrix 
% 
% Copyright (C) 2009, 2010 Carlo de Falco
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

function mat = op_div_v_q (spv, spq, msh)
  
  mat = spalloc(spq.ndof, spv.ndof, 1);
  
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      mat_loc = zeros (spq.nsh(iel), spv.nsh(iel));
      for idof = 1:spq.nsh(iel)
        ishp = spq.shape_functions(:, idof, iel);
        for jdof = 1:spv.nsh(iel) 
          jshd = spv.shape_function_divs(:, jdof, iel);
          %for inode = 1:msh.nqn
          mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
            sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* ishp .* jshd);
          %end  
        end
      end
      mat(spq.connectivity(:, iel), spv.connectivity(:, iel)) = ...
        mat(spq.connectivity(:, iel), spv.connectivity(:, iel)) + mat_loc;
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_div_v_q: singular map in element number %d', iel)
    end
  end

end
