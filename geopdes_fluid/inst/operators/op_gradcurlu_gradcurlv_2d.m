% OP_GRADCURLU_GRADCURLV: assemble the matrix A = [a(i,j)], a(i,j) = (mu grad (curl u_j), grad (curl v_i)).
%
%   mat = op_gradcurlu_gradcurlv (spu, spv, msh, mu);
%
% INPUT:
%    
%   spu:     structure representing the space of trial functions (see sp_bspline_2d_phys)
%   spv:     structure representing the space of test functions  (see sp_bspline_2d_phys)
%   msh:     structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%   mu:      coefficient
%
% OUTPUT:
%
%   mat: assembled matrix
% 
% Copyright (C) 2010 Rafael Vazquez, Carlo de Falco
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

function mat = op_gradcurlu_gradcurlv_2d (spu, spv, msh, coeff)
  
  mat = spalloc (spv.ndof, spu.ndof, 1);
  ndir = size (msh.quad_nodes, 1);
  
  gcu(1,1,:,:,:) =  spu.shape_function_hessians(2,1,:,:,:);
  gcu(1,2,:,:,:) =  spu.shape_function_hessians(2,2,:,:,:);
  gcu(2,1,:,:,:) = -spu.shape_function_hessians(1,1,:,:,:);
  gcu(2,2,:,:,:) = -spu.shape_function_hessians(1,2,:,:,:);

  gcv(1,1,:,:,:) =  spv.shape_function_hessians(2,1,:,:,:);
  gcv(1,2,:,:,:) =  spv.shape_function_hessians(2,2,:,:,:);
  gcv(2,1,:,:,:) = -spv.shape_function_hessians(1,1,:,:,:);
  gcv(2,2,:,:,:) = -spv.shape_function_hessians(1,2,:,:,:);

  for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
      mat_loc = zeros (spv.nsh(iel), spu.nsh(iel));
      for idof = 1:spv.nsh(iel)
        ishgc = reshape(gcv(:,:,:,idof,iel), ndir*ndir, []);
        for jdof = 1:spu.nsh(iel) 
          jshgc = reshape(gcu(:,:,:,jdof,iel), ndir*ndir, []);
          %for inode = 1:msh.nqn
          mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
             sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* ...
                  sum (ishgc .* jshgc, 1).' .* coeff(:,iel));
          %end  
        end
      end
      mat(spv.connectivity(1:spv.nsh(iel), iel), spu.connectivity(1:spu.nsh(iel), iel)) = ...
        mat(spv.connectivity(1:spv.nsh(iel), iel), spu.connectivity(1:spu.nsh(iel), iel)) + mat_loc;
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradcurlu_gradcurlv: singular map in element number %d', iel)
    end
  end

end
