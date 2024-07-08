% SP_GET_VERTEX_NEIGHBORS: Compute the indices of elements in the support of vertex functions, separated by patches.
%
% cell_indices = sp_get_vertex_neighbors (space, msh, vertex_index, patch_indices)
%
% INPUT:
%    space:   object defining the space of discrete functions (see sp_multipatch_C1)
%    msh:     object defining the domain partition and the quadrature rule (see msh_multipatch)
%    vertex_index: indices of the vertex for which to compute the cells.
%    patch_indices: local indices (w.r.t. vertex) for the patches on which to
%                    compute the cells.
%
% OUTPUT:
%    cell_indices: indices of the functions that interact with the given ones.
%
% Copyright (C) 2021, 2022 Rafael Vazquez
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

function cell_indices = sp_get_vertex_neighbors (space, msh, vertex_index, patch_indices)

if (nargin < 4)
  patch_indices = 1:space.vertices(vertex_index).valence_p;
end
cell_indices = cell (numel(patch_indices), 1);

Nelem = cumsum ([0 msh.nel_per_patch]);
for iptc = 1:numel(patch_indices)
  patch = space.vertices(vertex_index).patches(patch_indices(iptc));
  Cpatch_cols = sp_get_functions_on_patch (space, patch);
  [~,~,inds] = intersect (space.dofs_on_vertex{vertex_index}, Cpatch_cols);
  if (~isempty (inds))
    Cpatch = sp_compute_Cpatch (space, patch);
    [bsp_indices, ~] = find (Cpatch(:,inds));
    [aux_cell_indices, ~] = sp_get_cells (space.sp_patch{patch}, msh.msh_patch{patch}, bsp_indices(:)');
    cell_indices{iptc} = Nelem(patch) + aux_cell_indices(:).';
  else
    cell_indices{iptc} = [];
  end
end

end