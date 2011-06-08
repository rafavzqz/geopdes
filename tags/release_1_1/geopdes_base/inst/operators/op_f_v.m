% OP_F_V: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i).
%
%   mat = op_f_v (spv, msh, coeff);
%
% INPUT:
%     
%   spv:   structure representing the function space (see sp_bspline_2d_phys)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%   coeff: source function evaluated at the quadrature points
%
% OUTPUT:
%
%   mat: assembled right-hand side
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



function mat = op_f_v (spv, msh, coeff)
  
 mat   = zeros(spv.ndof, 1);

 shpv  = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel);
 coeff = reshape (coeff, spv.ncomp, msh.nqn, msh.nel);
 for iel = 1:msh.nel
   if (all (msh.jacdet(:,iel)))
     mat_loc = zeros (spv.nsh(iel), 1);
     for idof = 1:spv.nsh(iel)
       ishp = reshape(shpv(:,:,idof,iel), spv.ncomp, []);
       %        for inode = 1:msh.nqn
       mat_loc(idof) = mat_loc(idof) + ...
           sum (msh.jacdet(:, iel) .* msh.quad_weights(:, iel) .* ...
                sum(ishp .* coeff(:, :, iel), 1).');
       %        end  
     end
     mat(spv.connectivity(1:spv.nsh(iel), iel)) = mat(spv.connectivity(1:spv.nsh(iel), iel)) + mat_loc; 
   else
     warning ('geopdes:jacdet_zero_at_quad_node', 'op_f_v: singular map in element number %d', iel)
   end
 end
 
end
