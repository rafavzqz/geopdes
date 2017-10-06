% MSH_RESTRICT_TO_CELLS: restrict a msh struct to a subset of its elements
%  in elem_list.
%
% This function is useful for the current version of hierarchical splines. 
% The function allows to avoid recomputing the parametrization, which is
%  the most expensive part in the current version of hierarchical splines.
%
%     msh_elems = msh_restrict_to_cells (msh, elem_list)
%
% INPUTS:
%
%    msh:       mesh struct, computed with msh_cartesian/msh_evaluate_col(element_list)
%    elem_list: numbering of the elements, a subset of msh.elem_list.
%
% OUTPUT:
%
%     msh_elems: structure containing the quadrature rule in the given elements of the physical domain, which contains the following fields
%
%     FIELD_NAME         (SIZE)                  DESCRIPTION
%     ndim               (scalar)                dimension of the parametric space
%     rdim               (scalar)                dimension of the physical space
%     nel                (scalar)                number of elements in the sublist
%     elem_list          (1 x nel)               numbering of the elements in the sublist
%     nqn                (scalar)                number of quadrature points per element (must be the same for every patch)
%     nqn_dir            (1 x ndim)              number of quadrature points in each direction (must be the same for every patch)
%     quad_weights, geo_map, geo_map_jac, deo_map_der2, jacdet, element_size (see msh_evaluate_col for details)
%
% Copyright (C) 2017 Rafael Vazquez
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

function msh_col = msh_restrict_to_cells (msh, elems)

  msh_col.ndim = msh.ndim;
  msh_col.rdim = msh.rdim;

  [elem_list, ia, global_elem_list] = intersect (elems, msh.elem_list);
  if (numel (ia) ~= numel (elems))
    warning ('There are elements that are not in the original list')
  end
  msh_col.nel_dir = msh.nel_dir;
  msh_col.nel  = numel (elems);
  
  msh_col.nqn_dir = msh.nqn_dir;
  msh_col.nqn = msh.nqn;
  msh_col.elem_list = elem_list(:)';

  msh_col.quad_weights = msh.quad_weights(:,global_elem_list);
  msh_col.geo_map      = msh.geo_map(:,:,global_elem_list);
  msh_col.geo_map_jac  = msh.geo_map_jac(:,:,:,global_elem_list);
  msh_col.geo_map_der2 = msh.geo_map_der2(:,:,:,:,global_elem_list);
  msh_col.jacdet       = msh.jacdet(:,global_elem_list);
  msh_col.element_size = msh.element_size(:,global_elem_list);
  
end