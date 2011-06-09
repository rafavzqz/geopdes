% OP_F_V: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i).
%
%   rhs = op_f_v (spv, msh, coeff);
%
% INPUT:
%     
%   spv:   structure representing the function space (see sp_bspline_2d_phys)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%   coeff: source function evaluated at the quadrature points
%
% OUTPUT:
%
%   rhs: assembled right-hand side
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



function rhs = op_f_v2 (spv, msh, coeff, element_list)
  
 coeff = reshape (coeff, spv.ncomp, msh.nqn, msh.nel);

 if (nargin == 3)
   element_list = 1:msh.nel;
 end
 nel = numel (element_list);

 rhs   = zeros (spv.ndof, 1);
 shpv  = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, nel);

 for iel = 1:nel
   iel_glob = element_list (iel);
   if (all (msh.jacdet(:,iel_glob)))
     rhs_loc = zeros (spv.nsh(iel), 1);

     jacdet_weights = msh.jacdet(:, iel_glob) .* msh.quad_weights(:, iel_glob);
     jacdet_weights = repmat (jacdet_weights, [1, spv.nsh(iel)]);

     coeff_iel = repmat (coeff(:, :, iel_glob), [1, 1, spv.nsh(iel)]);

     shpv_iel = reshape (shpv(:, :, 1:spv.nsh(iel), iel), spv.ncomp, msh.nqn, spv.nsh(iel));

     val1 = reshape (sum (shpv_iel .* coeff_iel, 1), [msh.nqn, spv.nsh(iel)]);
     rhs_loc = sum (jacdet_weights .* val1, 1).';

     rhs(spv.connectivity(1:spv.nsh(iel), iel)) = rhs(spv.connectivity(1:spv.nsh(iel), iel)) + rhs_loc; 
   else
     warning ('geopdes:jacdet_zero_at_quad_node', 'op_f_v: singular map in element number %d', iel)
   end
 end
 
end
