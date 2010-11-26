% OP_GRADU_GRADV: assemble the stiffness matrix A = [a(i,j)], a(i,j) = (epsilon grad u_j, grad v_i).
%
%   mat = op_gradu_gradv (spu, spv, msh, epsilon);
%
% INPUT:
%    
%   spu:     structure representing the space of trial functions (see sp_bspline_2d_phys)
%   spv:     structure representing the space of test functions  (see sp_bspline_2d_phys)
%   msh:     structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%   epsilon: diffusion coefficient
%
% OUTPUT:
%
%   mat: assembled stiffness matrix
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

function mat = op_gradu_gradv (spu, spv, msh, coeff)
  
  mat = spalloc (spv.ndof, spu.ndof, 1);
  
  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);

  ndir = size (gradu, 2);

  for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
      mat_loc = zeros (spv.nsh(iel), spu.nsh(iel));
      for idof = 1:spv.nsh(iel)
        ishg = reshape(gradv(:,:,:,idof,iel),spv.ncomp * ndir, []);
        for jdof = 1:spu.nsh(iel) 
          jshg = reshape(gradu(:,:,:,jdof,iel),spu.ncomp * ndir, []);
          %for inode = 1:msh.nqn
          mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
             sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* ...
                  sum (ishg .* jshg, 1).' .* coeff(:,iel));
          %end  
        end
      end
      mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
        mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) + mat_loc;
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradu_gradv: singular map in element number %d', iel)
    end
  end

end
