% MSH_EVALUATE_ELEMENT_LIST: evaluate the parameterization in a given list of elements.
%
%     msh_elems = msh_evaluate_element_list (msh, elements)
%
% INPUTS:
%
%    msh:          mesh object (see msh_cartesian)
%    element_list: numbering of the elements where the evaluations are performed.
%
% OUTPUT:
%
%     msh_elems: structure containing the quadrature rule in the given elements of the physical domain, which contains the following fields
%
%     FIELD_NAME         (SIZE)                  DESCRIPTION
%     npatch             (scalar)                number of patches
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
% The function only works if the number of quadrature points is the same for all the patches and all directions.
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

function msh_col = msh_evaluate_element_list (msh, elem_list, varargin)

  elem_list = elem_list(:)';

  msh_col.npatch = msh.npatch;
  msh_col.ndim = msh.ndim;
  msh_col.rdim = msh.rdim;

  msh_col.nel  = numel (elem_list);
  msh_col.elem_list = elem_list;

  if (isempty (elem_list)), return, end
  
  fields = {'quad_weights', 'geo_map', 'geo_map_jac', 'geo_map_der2', 'jacdet', 'element_size'};
  cat_position = [2 3 4 5 2 2];
  for ii = 1:numel (fields)
    msh_col.(fields{ii}) = [];
  end
  
  Nelem = cumsum ([0, msh.nel_per_patch]);
  
  indices = cell (1, msh.npatch);
  for iptc = 1:msh.npatch
    [~,indices{iptc},~] = intersect ((Nelem(iptc)+1):Nelem(iptc+1), elem_list);
  end
  active_patches = find (~cellfun (@isempty, indices));

  
  for iptc = active_patches
    msh_patch = msh_evaluate_element_list (msh.msh_patch{iptc}, indices{iptc});
    
    if (msh_patch.nqn_dir ~= msh.msh_patch{active_patches(1)}.nqn_dir)
      error ('msh_evaluate_element_list: for multipatch geometries, the number of quadrature points should be the same on each patch')
    end

    for ii = 1:numel (fields)
      if (isfield (msh_patch, fields{ii}))
        msh_col.(fields{ii}) = cat (cat_position(ii), msh_col.(fields{ii}), msh_patch.(fields{ii}));
      end
    end
  end
  
  msh_col.nel_per_patch = cellfun (@numel, indices);
  msh_col.elem_list_of_patch = indices;
  msh_col.nel_dir_of_patch = cell (1, msh.npatch);
  for iptc = 1:msh.npatch
    msh_col.nel_dir_of_patch{iptc} = msh.msh_patch{iptc}.nel_dir;
  end
  msh_col.nqn = msh.msh_patch{active_patches(1)}.nqn;
  msh_col.nqn_dir = msh.msh_patch{active_patches(1)}.nqn_dir;

end