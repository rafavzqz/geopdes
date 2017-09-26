% SP_GET_BASIS_FUNCTIONS: Compute the indices of tensor-product B-splines acting on a list of cells.
%
% [fun_indices, indices_per_cell] = sp_get_basis_functions (space, msh, indices)
%
% INPUT:
%    space:   object defining the space of discrete functions (see sp_vector)
%    msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%    indices: indices of the cells.
%
% OUTPUT:
%    fun_indices: indices of the basis functions acting on the cells.
%    indices_per_cell: cell-array with indices of the basis functions on each cell.
%
% Copyright (C) 2015, 2016, 2017 Rafael Vazquez
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

function [function_indices, indices_per_cell] = sp_get_basis_functions (space, msh, cell_indices)

% space_aux = sp_precompute_param (space, msh, 'value', false);
% function_indices = space_aux.connectivity (:,cell_indices);
% function_indices = unique (function_indices(:));

indices = cell (space.ncomp_param, 1);
aux_per_cell = cell (space.ncomp_param, 1);
indices_per_cell = cell (numel (cell_indices), 1);
for icomp = 1:space.ncomp_param
  [indices_comp, aux_per_cell{icomp}] = sp_get_basis_functions (space.scalar_spaces{icomp}, msh, cell_indices);
  indices{icomp} = space.cumsum_ndof(icomp) + indices_comp;
  indices{icomp} = indices{icomp}(:);
end
function_indices = vertcat (indices{:});

if (nargout == 2)
  for icomp = 1:space.ncomp_param
    aux_per_cell{icomp} = cellfun (@(x) x + space.cumsum_ndof(icomp), aux_per_cell{icomp}, 'UniformOutput', false);
  end
  for iel = 1:numel(cell_indices)
    for icomp = 1:space.ncomp_param
      indices_per_cell{iel} = union (indices_per_cell{iel}, aux_per_cell{icomp}{iel});
    end
  end
end

end
