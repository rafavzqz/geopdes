% OP_F_VXN_3D: assemble the vector r = [r(i)], with  r(i) = (f, v_i x n), with n the normal exterior vector.
%
%   rhs = op_f_vxn_3d (spv, msh, coeff);
%
% INPUT:
%     
%   spv:   structure representing the function space (see sp_vector_3d_curl_transform/sp_eval_boundary_side)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_3d/msh_eval_boundary_side)
%   coeff: source function evaluated at the quadrature points
%
% OUTPUT:
%
%   rhs: assembled right-hand side
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



function rhs = op_f_vxn_3d (spv, msh, coeff)
  
  rhs = zeros(spv.ndof, 1);
  for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
      jacdet_weights = msh.jacdet(:, iel) .* msh.quad_weights(:, iel);

      rhs_loc = zeros (spv.nsh(iel), 1);
      for idof = 1:spv.nsh(iel)
        ishp = squeeze(spv.shape_functions(:,:,idof,iel));
        ishp_x_n = [ishp(2,:) .* msh.normal(3,:,iel) - ...
                       ishp(3,:) .* msh.normal(2,:,iel); ...
                    ishp(3,:) .* msh.normal(1,:,iel) - ...
                       ishp(1,:) .* msh.normal(3,:,iel); ...
                    ishp(1,:) .* msh.normal(2,:,iel) - ...
                       ishp(2,:) .* msh.normal(1,:,iel)];

          rhs_loc(idof) = rhs_loc(idof) +  ...
            sum (jacdet_weights .* sum(ishp_x_n .* coeff(:, :, iel), 1).');
      end

      rhs(spv.connectivity(1:spv.nsh(iel), iel)) = rhs(spv.connectivity(1:spv.nsh(iel), iel)) + rhs_loc;
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_f_vxn_3d: singular map in element number %d', iel)
    end
  end

end

