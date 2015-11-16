% SP_GET_CELLS: Compute the indices of the cells within the support of a list of tensor-product B-spline function.
%
% [cell_indices, indices_per_function] = sp_get_cells (space, msh, indices)
%
% INPUT:
%    space:   object defining the space of discrete functions (see sp_scalar)
%    msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%    indices: indices of the functions.
%
% OUTPUT:
%    cell_indices: indices of the cells within the support of the basis functions.
%    indices_per_function: indices of the cells within the support of each basis function.
%
% Copyright (C) 2015 Rafael Vazquez
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

function [cell_indices, indices_per_function] = sp_get_cells (space, msh, fun_indices)

space = sp_precompute_param (space, msh, 'value', false);

conn_indices = arrayfun (@(x) find (space.connectivity == x), fun_indices, 'UniformOutput', false);
[~, indices_per_function] = cellfun (@(x) ind2sub ([space.nsh_max, msh.nel], x), conn_indices, 'UniformOutput', false);
cell_indices = unique (vertcat (indices_per_function{:}));

% % Old version, to be used in case of unexpected memory problems
% subindices = cell (msh.ndim, 1);
% [subindices{:}] = ind2sub ([space.ndof_dir, 1], fun_indices); % The extra one makes it work in any dimension
% 
% indices_per_function = cell (numel (fun_indices), 1);
% for ifun = 1:numel (fun_indices)
%   cells = cell (msh.ndim, 1);
%   for idim = 1:msh.ndim
%     cells_1d{idim} = space.sp_univ(idim).supp{subindices{idim}(ifun)};
%   end
%   [cells{:}] = ndgrid (cells_1d{:});
%   indices_per_function{ifun} = sub2ind ([msh.nel_dir, 1], cells{:});
% end
% 
% indices_per_function = cellfun(@(x) x(:), indices_per_function, 'UniformOutput', false);
% cell_indices = unique (vertcat (indices_per_function{:}));

end
