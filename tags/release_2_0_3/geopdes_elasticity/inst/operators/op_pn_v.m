% OP_PN_V: assemble the right-hand side vector R = [r(i)], with  r(i) = (p n, v_i), where n is the normal vector.
%
%   mat = op_pn_v (spv, msh, press);
%
% INPUT:
%
%   spv:   structure representing the function space (see sp_bspline_2d/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_2d/msh_evaluate_col)
%   press: function evaluated at the quadrature points
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



function mat = op_pn_v (spv, msh, press)
  
 press = repmat (reshape (press, 1, msh.nqn, msh.nel), [spv.ncomp, 1, 1]);
 press = press .* msh.normal;
 mat   = op_f_v (spv, msh, press);

end
