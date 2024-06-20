% OP_F_curlV_3d: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, curl v_i).
%
%   rhs = op_f_curlv_3d (spv, msh, coeff);
%
% INPUT:
%
%   spv:   structure representing the function space (see sp_vector/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   coeff: source function evaluated at the quadrature points
%
% OUTPUT:
%
%   rhs: assembled right-hand side
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2017, 2019, 2023 Rafael Vazquez
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

function rhs = op_f_curlv_3d (spv, msh, coeff)
  
 coeff = reshape (coeff, spv.ncomp, [], msh.nqn, msh.nel);

 rhs   = zeros (spv.ndof, 1);
 curlv  = reshape (spv.shape_function_curls, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel);

 jacdet_weights = msh.jacdet .* msh.quad_weights;

 for iel = 1:msh.nel
   if (all (msh.jacdet(:,iel)))
     coeff_iel = reshape (coeff(:,:,:,iel), spv.ncomp, msh.nqn, 1, 1);
     jacdet_iel = reshape (jacdet_weights(:, iel), [1, msh.nqn, 1, 1]);
     coeff_times_jw = bsxfun (@times, jacdet_iel, coeff_iel);

     curlv_iel = reshape (curlv(:,:,:,iel), spv.ncomp, msh.nqn, spv.nsh_max, 1);

     aux_val = bsxfun (@times, coeff_times_jw, curlv_iel);
     rhs_loc = sum (sum (aux_val, 1), 2);

     indices = find (spv.connectivity(:,iel));
     rhs_loc = rhs_loc(indices); conn_iel = spv.connectivity(indices,iel);
     rhs(conn_iel) = rhs(conn_iel) + rhs_loc(:); 
   else
     warning ('geopdes:jacdet_zero_at_quad_node', 'op_f_curlv_3d: singular map in element number %d', iel)
   end
 end
 
end
