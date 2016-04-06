% MSH_RESTRICT_TO_PATCHES: restrict a struct, computed for a given list of
%   elements, to the elements that only belong to a certain list of patches.
%
% This function is useful for the current version of hierarchical splines, 
%  in particular for boundary conditions, since only some patches need to
%  be computed.
% The function allows to avoid recomputing the parametrization, which is
%  the most expensive part in the current version of hierarchical splines.
%
%     msh_elems = msh_restrict_to_patches (msh, patch_list)
%
% INPUTS:
%
%    msh:        mesh struct, computed with msh_multipatch/msh_evaluate_element_list
%    patch_list: numbering of the elements where the evaluations are performed.
%
% OUTPUT:
%
%     msh_elems: structure containing the quadrature rule in the given elements of the physical domain, which contains the following fields
%
%     FIELD_NAME         (SIZE)                  DESCRIPTION
%     npatch             (scalar)                total number of patches
%     ndim               (scalar)                dimension of the parametric space
%     rdim               (scalar)                dimension of the physical space
%     nel                (scalar)                number of elements in the list
%     elem_list          (1 x nel)               numbering of the elements in the list
%     nqn                (scalar)                number of quadrature points per element (must be the same for every patch)
%     nqn_dir            (1 x ndim)              number of quadrature points in each direction (must be the same for every patch)
%     nel_per_patch      (1 x npatch)            number of selected elements on each patch
%     elem_list_of_patch (1 x npatch cell-array) selected elements on the patch, with local numbering
%     nel_dir_of_patch   (1 x npatch cell-array) the total number of elements in each direction, for each patch
%     quad_weights, geo_map, geo_map_jac, deo_map_der2, jacdet, element_size (see msh_evaluate_col for details)
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

function msh_col = msh_restrict_to_patches (msh, patches)

  msh_col.npatch = msh.npatch;
  msh_col.ndim = msh.ndim;
  msh_col.rdim = msh.rdim;

  msh_col.nel  = [];
  
  msh_col.nel_per_patch = zeros (1, msh.npatch);
  msh_col.nel_dir_of_patch = cell (1, msh.npatch);
  msh_col.nel_per_patch(patches) = msh.nel_per_patch(patches);
  msh_col.nel_dir_of_patch = msh.nel_dir_of_patch;
  msh_col.nel = sum (msh_col.nel_per_patch);

  msh_col.elem_list = [];
  
  nonactive_patches = setdiff (1:msh.npatch, patches);
  msh_col.elem_list_of_patch = msh.elem_list_of_patch;
  for iptc = nonactive_patches
    msh_col.elem_list_of_patch{iptc} = [];
  end
  msh_col.nqn = msh.nqn;
  msh_col.nqn_dir = msh.nqn_dir;
  
  Nelem = cumsum ([0, msh.nel_per_patch]);
  global_elem_list = [];
  for iptc = patches
    global_elem_list = union (global_elem_list, Nelem(iptc)+1:Nelem(iptc+1));
  end
  msh_col.elem_list = msh.elem_list(global_elem_list);

  msh_col.quad_weights = msh.quad_weights(:,global_elem_list);
  msh_col.geo_map      = msh.geo_map(:,:,global_elem_list);
  msh_col.geo_map_jac  = msh.geo_map_jac(:,:,:,global_elem_list);
  msh_col.geo_map_der2 = msh.geo_map_der2(:,:,:,:,global_elem_list);
  msh_col.jacdet       = msh.jacdet(:,global_elem_list);
  msh_col.element_size = msh.element_size(:,global_elem_list);
  
end