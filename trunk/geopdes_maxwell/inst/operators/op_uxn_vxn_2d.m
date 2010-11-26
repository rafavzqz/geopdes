% OP_UXN_VXN_2D: assemble the matrix M = [m(i,j)], m(i,j) = (mu u_j x n, v_i x n), with n the exterior normal vector.
%
%   mat = op_uxn_vxn_2d (spu, spv, msh, coeff);
%
% INPUT:
%   
%  spu:   structure representing the space of trial functions (see sp_bsp_hcurl_2d_phys)
%  spv:   structure representing the space of test functions  (see sp_bsp_hcurl_2d_phys)
%  msh:   structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%  coeff: physical parameter
%
% OUTPUT:
%
%  mat: assembled matrix
% 
% Copyright (C) 2009, 2010 Carlo de Falco, Rafael Vazquez
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

function mat = op_uxn_vxn_2d (spu, spv, msh, coeff)
  
  mat = spalloc(spv.ndof, spu.ndof, 1);
  for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
      mat_loc = zeros (spv.nsh(iel), spu.nsh(iel));
      for idof = 1:spv.nsh(iel)
        ishp = squeeze(spv.shape_functions(:,:,idof,iel));
        ishp_x_n = (ishp(1,:) .* msh.normal(2,:,iel) - ...
                    ishp(2,:) .* msh.normal(1,:,iel))';
        for jdof = 1:spu.nsh(iel)
          jshp = squeeze(spu.shape_functions(:,:,jdof,iel));
          jshp_x_n = (jshp(1,:) .* msh.normal(2,:,iel) - ...
                      jshp(2,:) .* msh.normal(1,:,iel))';
          %for inode = 1:msh.nqn
            mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
              sum (msh.jacdet(:, iel) .* msh.quad_weights(:, iel) .* ...
              ishp_x_n .* jshp_x_n .* coeff(:, iel));
          %end  
        end
      end
      mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
        mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) + mat_loc;
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_uxn_vxn_2d: singular map in element number %d', iel)
    end
  end

end

