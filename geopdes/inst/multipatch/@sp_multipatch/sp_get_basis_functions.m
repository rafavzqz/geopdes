% SP_GET_BASIS_FUNCTIONS: Compute the indices of tensor-product B-splines acting on a list of cells.
%
% fun_indices = sp_get_basis_functions (space, msh, indices)
%
% INPUT:
%    space:   object defining the space of discrete functions (see sp_multipatch)
%    msh:     object defining the domain partition and the quadrature rule (see msh_multipatch)
%    indices: indices of the cells.
%
% OUTPUT:
%    fun_indices: indices of the basis functions acting on the cells.
%
% Copyright (C) 2015, 2016 Rafael Vazquez
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

function function_indices = sp_get_basis_functions (space, msh, cell_indices)

function_indices = [];
Nelem = cumsum ([0 msh.nel_per_patch]);
for iptc = 1:space.npatch
  [~,local_cell_indices,~] = intersect ((Nelem(iptc)+1):Nelem(iptc+1), cell_indices);
  if (~isempty (local_cell_indices))
    aux_indices = sp_get_basis_functions (space.sp_patch{iptc}, msh.msh_patch{iptc}, local_cell_indices);
    function_indices = union (function_indices, space.gnum{iptc}(aux_indices));
  end
end

end
