% MSH_CELLS_NEAR_VERTEX: Identify the elements adjacent to given vertices.
%
%     [all_adjacent_elemens, elems_near_vertex] = msh_cells_near_vertex (msh, vertices)
%
% INPUT:
%
%    msh:      multipatch mesh object (see msh_multipatch)
%    vertices: struct containing the vertex information (see vertices_struct)
%
% OUTPUT:
%
%    all_adjacent_elements: elements adjacents to any of the input vertices.
%    elems_near_vertex: for each vertex in the input, list of elements
%        adjacent to the vertex (cell-array, to allow different valences).
%
% Copyright (C) 2021 Rafael Vazquez
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

function [all_adj_elems, elems_adj_to_vertex] = msh_cells_near_vertex (msh, vertices)

nvert = numel (vertices);
elems_adj_to_vertex = cell (nvert, 1);
shifting_index = cumsum ([0 msh.nel_per_patch]);

all_adj_elems = [];
for ivert = 1:nvert
  elems = [];
  patches = vertices(ivert).patches;
  patch_ornt = vertices(ivert).patch_reorientation;
  for iptc = 1:numel(patches)
    patch = patches(iptc);
    nel_dir = msh.msh_patch{patch}.nel_dir;
    indx = 1; indy = 1;
    if (patch_ornt(iptc,1))
      indx = nel_dir(1);
    end
    if (patch_ornt(iptc,2))
      indy = nel_dir(2);
    end
    new_elem = shifting_index(patch) + sub2ind (nel_dir, indx, indy);
    elems = [elems, new_elem];
  end
  elems_adj_to_vertex{ivert} = elems;
  all_adj_elems = union (all_adj_elems, elems);
end

end