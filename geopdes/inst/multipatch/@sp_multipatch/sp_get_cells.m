% SP_GET_CELLS: Compute the indices of the cells within the support of a list of tensor-product B-spline function.
%
% [cell_indices, indices_per_function] = sp_get_cells (space, msh, indices)
%
% INPUT:
%    space:   object defining the space of discrete functions (see sp_multipatch)
%    msh:     object defining the domain partition and the quadrature rule (see msh_multipatch)
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

fun_indices = transpose (fun_indices(:));

indices_per_function = cell (numel (fun_indices), 1);
cell_indices = [];
Nelem = cumsum ([0 msh.nel_per_patch]);
for iptc = 1:space.npatch
  sp_patch = sp_precompute_param (space.sp_patch{iptc}, msh.msh_patch{iptc}, 'value', false);
  sp_patch.connectivity = space.gnum{iptc}(sp_patch.connectivity);
  
  conn_indices = arrayfun (@(x) find (sp_patch.connectivity == x), fun_indices, 'UniformOutput', false);
  [~, ind_per_fun] = cellfun (@(x) ind2sub ([sp_patch.nsh_max, msh.msh_patch{iptc}.nel], x), conn_indices, 'UniformOutput', false);
  
  ind_per_fun = cellfun (@(x) x + Nelem(iptc), ind_per_fun, 'UniformOutput', false);
  cell_indices = union (cell_indices, vertcat (ind_per_fun{:}));

  if (nargout == 2)
    [~,local_funs,~] = intersect (fun_indices, space.gnum{iptc});
    for ifun = local_funs
      indices_per_function{ifun} = union (indices_per_function{ifun}, ind_per_fun{ifun});
    end
  end
end
