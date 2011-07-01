% SP_EVAL_BOUNDARY_SIDE: Construct the space structure of one side of the boundary.
%
%     sp = sp_eval_boundary_side (sp, msh, iside)
%
% INPUTS:
%     
%
% OUTPUT:
%
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

function sp_bnd = sp_eval_boundary_side (sp, msh_side)

  iside = msh_side.side_number;
  ind = mod (floor ((iside+1)/2), 2) + 1;

  bnodes = reshape (squeeze (msh_side.quad_nodes(ind,:,:)), msh_side.nqn, []);
  sp_bnd = sp_bspline_1d_param (sp.knots{ind}, sp.degree(ind), bnodes, ...
                                'gradient', false);

  sp_bnd.dofs = sp.boundary(iside).dofs;

end
