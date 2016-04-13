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

indices_per_comp = cell (numel (space.ncomp_param), 1);
for icomp = 1:space.ncomp_param
  [~,indices,~] = intersect (space.cumsum_ndof(icomp)+1:space.cumsum_ndof(icomp+1), fun_indices);
  
  subindices = cell (msh.ndim, 1);
  [subindices{:}] = ind2sub ([space.scalar_spaces{icomp}.ndof_dir, 1], indices); % The extra one makes it work in any dimension

  indices_per_function = cell (numel (indices), 1);
  for ifun = 1:numel(indices)
    funs = cell (msh.ndim, 1);
    fun_1d = cell (msh.ndim, space.ncomp_param);
    elem_1d = cell (msh.ndim, 1);
    indaux = cell (space.ncomp_param, 1);
    
    for idim = 1:msh.ndim
      [~,elem_1d{idim}] = find (space.scalar_spaces{icomp}.sp_univ(idim).connectivity == subindices{idim}(ifun));
      elem_1d{idim} = unique (elem_1d{idim});
    end
    
    for jcomp = 1:space.ncomp_param
      for idim = 1:msh.ndim
        fun_1d{idim,jcomp} = unique (space.scalar_spaces{jcomp}.sp_univ(idim).connectivity(:,elem_1d{idim}));
      end
      [funs{:}] = ndgrid (fun_1d{:,jcomp});
      indaux{jcomp} = space.cumsum_ndof(jcomp) + unique (sub2ind ([space.scalar_spaces{jcomp}.ndof_dir, 1], funs{:}));
    end
  
    indices_per_function{ifun} = vertcat (indaux{:});
  end
  indices_per_comp{icomp} = unique (vertcat (indices_per_function{:}));
end

neighbors_indices = unique (vertcat (indices_per_comp{:}));

end