% SP_GET_BOUNDARY_FUNCTIONS: indices of the degrees of freedom of the given boundaries.
%
%   dofs = sp_get_boundary_functions (space, msh, sides)
%
% INPUT:
%
%  space:  object representing the space of trial functions (see sp_multipatch)
%  msh:    object containing the domain partition and the quadrature rule (see msh_multipatch)
%  sides:  boundary sides from which we want to compute the indices
%
% OUTPUT:
%
%  dofs: global numbering of the boundary basis functions
%
% Copyright (C) 2019 Rafael Vazquez
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

function dofs = sp_get_boundary_functions (space, msh, refs)

  boundaries = msh.boundaries;
  Nbnd = cumsum ([0, boundaries.nsides]);
  bnd_dofs = [];
  for iref = refs
    iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    boundary_gnum = space.boundary.gnum;
    bnd_dofs = union (bnd_dofs, [boundary_gnum{iref_patch_list}]);
  end
  
  dofs = space.boundary.dofs(bnd_dofs);
  
end
