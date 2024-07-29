% SP_GET_CELLS: Compute the indices of the cells within the support of a list of B-spline functions.
%
% [cell_indices, indices_per_function] = sp_get_cells (space, msh, indices)
%
% INPUT:
%    space:   object defining the space of discrete functions (see sp_multipatch_C1)
%    msh:     object defining the domain partition and the quadrature rule (see msh_multipatch)
%    indices: indices of the functions.
%
% OUTPUT:
%    cell_indices: indices of the cells within the support of the basis functions.
%    indices_per_function: indices of the cells within the support of each basis function.
%
% Copyright (C) 2015, 2016, 2022 Rafael Vazquez
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

fun_indices = fun_indices(:).';

indices_per_function = cell (numel (fun_indices), 1);
cell_indices = [];

Nelem = cumsum ([0 msh.nel_per_patch]);
for iptc = 1:space.npatch
  Cpatch_cols = sp_get_functions_on_patch (space, iptc);
  [~,loc_f_inds,fun_indices_on_patch] = intersect (fun_indices, Cpatch_cols);
  if (~isempty (fun_indices_on_patch))
    Cpatch = sp_compute_Cpatch (space, iptc);
    [patch_indices,local_funs] = find (Cpatch(:,fun_indices_on_patch));
    patch_indices = patch_indices(:).';
    [aux_cell_indices, ind_per_fun] = sp_get_cells (space.sp_patch{iptc}, msh.msh_patch{iptc}, patch_indices);

    cell_indices = union (cell_indices, Nelem(iptc)+aux_cell_indices);

    if (nargout == 2)
      local_funs = local_funs(:).';
      for ifun = 1:numel(patch_indices)
        indices_per_function{loc_f_inds(local_funs(ifun))} = union (indices_per_function{loc_f_inds(local_funs(ifun))}, Nelem(iptc)+ind_per_fun{ifun});
      end
    end
  end
end

end