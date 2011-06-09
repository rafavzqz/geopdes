% OP_F_V_TP: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i), using the tensor product structure.
%
%   rhs = op_f_v_tp (spv, msh, coeff);
%
% INPUT:
%     
%   spv:   class representing the function space (see sp_bspline_2d)
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

function rhs = op_f_v_tp (space, msh, coeff)

  rhs = zeros (space.ndof, 1);
  for iel = 1:msh.nelu
    [sp, elem_list] = sp_evaluate_col (space, msh, iel, 'gradient', false);
    rhs = rhs + op_f_v (sp, msh, coeff, elem_list);
  end

end
