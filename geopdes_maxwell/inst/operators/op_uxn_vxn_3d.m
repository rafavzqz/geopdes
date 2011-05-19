% OP_UXN_VXN_3D: assemble the matrix M = [m(i,j)], m(i,j) = (mu u_j x n, v_i x n), with n the exterior normal vector.
%
%   mat = op_uxn_vxn_3d (spu, spv, msh, coeff);
%
% INPUT:
%   
%  spu:   structure representing the space of trial functions (see sp_bsp_hcurl_3d_phys)
%  spv:   structure representing the space of test functions  (see sp_bsp_hcurl_3d_phys)
%  msh:   structure containing the domain partition and the quadrature rule (see msh_push_forward_3d)
%  coeff: physical parameter
%
% OUTPUT:
%
%  mat: assembled matrix
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

function mat = op_uxn_vxn_3d (spu, spv, msh, coeff)
  
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      jacdet_weights = msh.jacdet(:, iel) .* ...
              msh.quad_weights(:, iel) .* coeff(:, iel);
      for idof = 1:spv.nsh(iel)
        ishp = squeeze (spv.shape_functions(:, :, idof, iel));
        ishp_x_n = [ishp(2, :) .* msh.normal(3, :, iel) - ...
                       ishp(3, :) .* msh.normal(2, :, iel); ...
                    ishp(3, :) .* msh.normal(1, :, iel) - ...
                       ishp(1, :) .* msh.normal(3, :, iel); ...
                    ishp(1, :) .* msh.normal(2, :, iel) - ...
                       ishp(2, :) .* msh.normal(1, :, iel)];
        for jdof = 1:spu.nsh(iel)
          ncounter = ncounter + 1;
          rows(ncounter) = spv.connectivity(idof, iel);
          cols(ncounter) = spu.connectivity(jdof, iel);

          jshp = squeeze (spu.shape_functions(:, :, jdof, iel));
          jshp_x_n = [jshp(2, :) .* msh.normal(3, :, iel) - ...
                         jshp(3, :) .* msh.normal(2, :, iel); ...
                      jshp(3, :) .* msh.normal(1, :, iel) - ...
                         jshp(1, :) .* msh.normal(3, :, iel); ...
                      jshp(1, :) .* msh.normal(2, :, iel) - ...
                         jshp(2, :) .* msh.normal(1, :, iel)];

          values(ncounter) = ...
                 sum (jacdet_weights .* sum(ishp_x_n .* jshp_x_n, 1)');
        end
      end
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_uxn_vxn_3d: singular map in element number %d', iel)
    end
  end

  mat = sparse (rows, cols, values, spv.ndof, spu.ndof);

end

