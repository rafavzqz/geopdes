% SP_GET_BASIS_FUNCTIONS: Compute the indices of tensor-product B-splines acting on a list of cells.
%
% fun_indices = sp_get_basis_functions (space, msh, indices)
%
% INPUT:
%    space:   object defining the space of discrete functions (see sp_vector)
%    msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%    indices: indices of the cells.
%
% OUTPUT:
%    fun_indices: indices of the basis functions acting on the cells.
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

function function_indices = sp_get_basis_functions (space, msh, cell_indices)

% space_aux = sp_precompute_param (space, msh, 'value', false);
% function_indices = space_aux.connectivity (:,cell_indices);
% function_indices = unique (function_indices(:));
% 
subindices = cell (msh.ndim, 1);
[subindices{:}] = ind2sub ([msh.nel_dir, 1], cell_indices); % The extra one makes it work in any dimension

indices = cell (numel (cell_indices), space.ncomp_param);

for iel = 1:numel (cell_indices)
  for icomp = 1:space.ncomp_param
    conn = cell (msh.ndim, 1);
    conn_1d = cell (msh.ndim, 1);
    for idim = 1:msh.ndim
      conn_1d{idim} = space.scalar_spaces{icomp}.sp_univ(idim).connectivity(:,subindices{idim}(iel));
    end
    [conn{:}] = ndgrid (conn_1d{:});

    indices{iel,icomp} = space.cumsum_ndof(icomp) + sub2ind ([space.scalar_spaces{icomp}.ndof_dir, 1], conn{:}); % The extra 1 makes things work in any dimension
  end
end
function_indices = [indices{:}];
function_indices = unique (function_indices(:));

end
