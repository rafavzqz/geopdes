% SP_GET_NEIGHBORS: Compute the indices of functions that share one element in the support of a given list of functions
%
% neighbors_indices = sp_get_neighbors (space, msh, indices)
%
% INPUT:
%    space:   object defining the space of discrete functions (see sp_vector)
%    msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%    indices: indices of the functions.
%
% OUTPUT:
%    neighbors_indices: indices of the functions that interact with the given ones.
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

function neighbors_indices = sp_get_neighbors (space, msh, fun_indices)

% space_aux = sp_precompute_param (space, msh, 'value', false);
% conn_indices = arrayfun (@(x) find (space_aux.connectivity == x), fun_indices, 'UniformOutput', false);
% [~, indices_per_function] = cellfun (@(x) ind2sub ([space_aux.nsh_max, msh.nel], x), conn_indices, 'UniformOutput', false);
% cell_indices = unique (vertcat (indices_per_function{:}));
% 
% neighbors_indices = unique (space_aux.connectivity (:,cell_indices));

cells = sp_get_cells (space, msh, fun_indices);
neighbors_indices = sp_get_basis_functions (space, msh, cells);

end