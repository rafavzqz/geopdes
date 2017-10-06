% SP_GET_BASIS_FUNCTIONS: Compute the indices of tensor-product B-splines acting on a list of cells.
%
% [fun_indices, indices_per_cell] = sp_get_basis_functions (space, msh, indices)
%
% INPUT:
%    space:   object defining the space of discrete functions (see sp_scalar)
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

% space = sp_precompute_param (space, msh, 'value', false);
% function_indices = space.connectivity (:,cell_indices);
% function_indices = unique (function_indices(:));

% Old version, to be used in case of unexpected memory problems
ndim = msh.ndim;
ndof_dir = space.ndof_dir;
sp_univ = space.sp_univ;

subindices = cell (ndim, 1);
[subindices{:}] = ind2sub ([msh.nel_dir, 1], cell_indices); % The extra one makes it work in any dimension

indices_per_cell = cell (numel (cell_indices), 1);
conn = cell (ndim, 1);
conn_1d = cell (ndim, 1);
for iel = 1:numel (cell_indices)
  for idim = 1:ndim
    conn_1d{idim} = sp_univ(idim).connectivity(:,subindices{idim}(iel));
  end
  [conn{:}] = ndgrid (conn_1d{:});

  indices_per_cell{iel} = sub2ind ([ndof_dir, 1], conn{:}); % The extra 1 makes things work in any dimension
end
function_indices = [indices_per_cell{:}];
function_indices = unique (function_indices(:));

if (nargout == 2)
  indices_per_cell = cellfun (@(x) x(:), indices_per_cell, 'UniformOutput', false);
end

end
