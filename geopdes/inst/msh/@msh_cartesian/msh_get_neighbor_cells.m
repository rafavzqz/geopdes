% MSH_GET_NEIGHBOR_CELLS: Compute the indices of the neighbor cells of a given cell
%
% [cell_indices, indices_per_face] = msh_get_neighbor_cells (msh, indices, [diagonal])
%
% INPUT:
%    msh:      object defining the domain partition (see msh_cartesian)
%    indices:  indices of the input cells.
%    diagonal: decide whether to include the neighbors in diagonal direction (false by default)
%
% OUTPUT:
%    cell_indices: indices of the neighbor elements in global numbering.
%    indices_per_cell: indices of the neighbor cells for each input cell,
%      in an array of size (numel(indices), 2*ndim) when diagonal is false,
%      and an array of size (numel(indices), 3^ndim-1) when diagonal is true
%
% Copyright (C) 2017, 2018 Rafael Vazquez
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

function [neig_indices, indices_per_cell] = msh_get_neighbor_cells (msh, cell_indices, diagonal)

if (nargin < 3)
  diagonal = false;
end

subindices = cell (msh.ndim, 1);
indices_per_cell = zeros (numel (cell_indices), 2*msh.ndim);
neig_indices = [];

[subindices{:}] = ind2sub ([msh.nel_dir, 1], cell_indices(:)');

if (~diagonal || msh.ndim == 1)
  for idim = 1:msh.ndim
    subinds = subindices;
    
% Indices in the negative direction
    subinds{idim} = subindices{idim} - 1;
    index_out = find (subinds{idim} == 0);
    index_in = setdiff (1:numel(cell_indices), index_out);
    subinds_tmp = cellfun (@(x) x(index_in), subinds, 'UniformOutput', false);

    cell_inds = sub2ind ([msh.nel_dir, 1], subinds_tmp{:});
    indices_per_cell(index_out, 2*idim-1) = 0;
    indices_per_cell(index_in, 2*idim-1) = cell_inds;
    neig_indices = union (neig_indices, cell_inds);

% Indices in the positive direction
    subinds{idim} = subindices{idim} + 1;
    index_out = find (subinds{idim} == msh.nel_dir(idim) + 1);
    index_in = setdiff (1:numel(cell_indices), index_out);
    subinds_tmp = cellfun (@(x) x(index_in), subinds, 'UniformOutput', false);

    cell_inds = sub2ind ([msh.nel_dir, 1], subinds_tmp{:});
    indices_per_cell(index_out,2*idim) = 0;
    indices_per_cell(index_in, 2*idim) = cell_inds;
    neig_indices = union (neig_indices, cell_inds);

  end
  
else
  indices_per_cell = zeros (numel(cell_indices), 3^msh.ndim-1);
  indices_in = cell (msh.ndim, 1);
  subs_cell = cell (msh.ndim, 1);
  all_neigs = cell (msh.ndim, 1);
  
  mid_index = (3^msh.ndim + 1)/2;
  for iel = 1:numel (cell_indices)
    index_in = true;
    for idim = 1:msh.ndim
      subs_idim = subindices{idim}(iel) + [-1 0 1];
      indices_in{idim} = subs_idim > 0 & subs_idim < msh.nel_dir(idim)+1;
      subs_cell{idim} = subs_idim(indices_in{idim});
      index_in = kron (indices_in{idim}, index_in);
    end
    [all_neigs{:}] = ndgrid (subs_cell{:});
    neighbors = sub2ind ([msh.nel_dir, 1], all_neigs{:});
    neighbors = setdiff (neighbors, cell_indices(iel), 'stable');
    index_in = [index_in(1:mid_index-1) index_in(mid_index+1:end)];

    indices_per_cell(iel,logical(index_in)) = neighbors;
    neig_indices = union (neig_indices, neighbors);
  end
  
end
if iscolumn(neig_indices)
    neig_indices = neig_indices';
end
end
