% SP_GET_FUNCTIONS_ON_PATCH: Compute the indices of non-vanishing basis functions on a patch.
%
% [Cpatch_cols, dofs_interior, dofs_edge, dofs_vertex] = sp_get_functions_on_patch (space, patch)
%
% INPUT:
%    space:   object defining the space of discrete functions (see sp_multipatch_C1)
%    patch:   index of the patch
%
% OUTPUT:
%    Cpatch_cols: indices of the C^1 basis functions that do not vanish on the patch.
%    dofs_interior: indices of the C^1 interior basis functions.
%    dofs_edge: indices of the C^1 edge basis functions.
%    dofs_vertex: indices of the C^1 vertex basis functions.
%
% Copyright (C) 2022 Cesare Bracco, Rafael Vazquez
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

function [Cpatch_cols, dofs_interior, dofs_edge, dofs_vertex] = sp_get_functions_on_patch (space, iptc_ind)

  numel_interior_dofs = space.ndof_interior_per_patch;
  ndof_per_interface = cellfun (@numel, space.dofs_on_edge);
  ndof_per_vertex = cellfun (@numel, space.dofs_on_vertex);

% Interior basis functions
  global_indices = sum (numel_interior_dofs(1:iptc_ind-1)) + (1:numel_interior_dofs(iptc_ind));
  dofs_interior = global_indices;

% Edge basis functions
  dofs_edge = [];
  for intrfc = 1:numel(space.interfaces)
    patches = [space.interfaces(intrfc).patch1 space.interfaces(intrfc).patch2];
    if (ismember (iptc_ind, patches))
      global_indices = space.ndof_interior + sum(ndof_per_interface(1:intrfc-1)) + (1:ndof_per_interface(intrfc));
      dofs_edge = union (dofs_edge, global_indices);
    end
  end

% Vertex basis functions
  dofs_vertex = [];
  for ivrt = 1:numel(space.vertices)
    if (ismember (iptc_ind, space.vertices(ivrt).patches))
      global_indices = space.ndof_interior + space.ndof_edges + sum(ndof_per_vertex(1:ivrt-1)) + (1:ndof_per_vertex(ivrt));
      dofs_vertex = union (dofs_vertex, global_indices);
    end
  end
  
  Cpatch_cols = union (union (dofs_interior, dofs_edge), dofs_vertex);

end
