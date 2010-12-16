% OP_V_GRADP: assemble the matrix B = [b(i,j)], b(i,j) = (epsilon grad p_i, v_j).
%
%   mat = op_v_gradp (spv, spp, msh, epsilon);
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

function mat = op_v_gradp (spv, spp, msh, coeff)

  mat = spalloc(spp.ndof, spv.ndof, 1);
  ndir = size (spp.shape_function_gradients, 1);
  for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
      mat_loc = zeros (spp.nsh(iel), spv.nsh(iel));
      for idof = 1:spp.nsh(iel)
        ishg = reshape(spp.shape_function_gradients(:,:,idof,iel),ndir,[]);
        for jdof = 1:spv.nsh(iel)
          jshg = reshape(spv.shape_functions(:,:,jdof,iel),ndir,[]);
          %for inode = 1:msh.nqn
            mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
             sum (msh.jacdet(:, iel) .* msh.quad_weights(:, iel) .* ...
             sum (ishg .* jshg, 1).' .* coeff(:, iel));
          %end  
        end
      end
      mat(spp.connectivity(1:spp.nsh(iel), iel), spv.connectivity(1:spv.nsh(iel), iel)) = ...
        mat(spp.connectivity(1:spp.nsh(iel), iel), spv.connectivity(1:spv.nsh(iel), iel)) + mat_loc;
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_v_gradp: singular map in element number %d', iel)
    end
  end

end

