% SP_EVALUATE_ELEMENT_LIST: compute the basis functions in a given list of elements.
%
%     sp = sp_evaluate_element_list (space, msh_elems, 'option1', value1, ...)
%
% INPUTS:
%     
%    space:     object defining the space of discrete functions (see sp_multipatch)
%    msh_elems: structure containing the information of quadrature or
%               visualization points, for a given list of elements and patches
%               (see msh_multipatch/msh_evaluate_element_list)
%   'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
%            gradient   |      false      |  compute shape_function_gradients
%            hessian    |      false      |  compute shape_function_hessians (only scalars)
%            laplacian  |      false      |  compute shape_function_laplacians (only scalars)
%            div        |      false      |  compute shape_function_divs (only vectors)
%            curl       |      false      |  compute shape_function_curls (only vectors)
%
% OUTPUT:
%
%    sp: struct representing the discrete function space, with the following fields:
%
%    FIELD_NAME      (SIZE)                    DESCRIPTION
%    npatch          (scalar)                  number of patches
%    ncomp           (scalar)                  number of components of the functions of the space
%    ndof            (scalar)                  total number of degrees of freedom
%    nsh             (1 x nel)                 number of non-vanishing functions on each element
%    nsh_max         (scalar)                  maximum number of nsh
%    connectivity    (nsh_max x nel)           global numbering of the non-vanishing functions on each element
%    shape_functions, shape_function_gradients... (see sp_evaluate_col for details)
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

function sp = sp_evaluate_element_list (space, msh, varargin)

  sp.npatch = space.npatch;
  sp.ncomp = space.ncomp;
  sp.ndof = space.ndof;

  if (isempty (msh.elem_list)), return, end

  sp.nsh_max = [];
  fields = {'nsh', 'connectivity', 'shape_functions', 'shape_function_gradients', 'shape_function_hessians', 'shape_function_laplacians'};
  cat_position = [2, 2, 3, 4, 5, 3];
  for ii = 1:numel(fields)
    sp.(fields{ii}) = [];
  end
  
  for iptc = 1:space.npatch
%     msh_patch = msh_evaluate_element_list (msh.msh_patch{iptc}, msh.elem_list_of_patch{iptc}); % Not working anymore
    msh_patch = msh_restrict_to_patch (msh, iptc);

    sp_patch = sp_evaluate_element_list (space.sp_patch{iptc}, msh_patch, varargin{:});
    sp_patch.connectivity = reshape (space.gnum{iptc}(sp_patch.connectivity), sp_patch.nsh_max, msh_patch.nel);

    for ii = 1:numel(fields)
      if (isfield (sp_patch, fields{ii}))
        sp.(fields{ii}) = cat (cat_position(ii), sp.(fields{ii}), sp_patch.(fields{ii}));
      end
    end
  end
  sp.nsh_max = max (sp.nsh);

  for ii = 1:numel(fields)
    if (isempty (sp.(fields{ii})))
      sp = rmfield (sp, fields{ii});
    end
  end
  
  if (~isempty (space.dofs_ornt))
    error ('You have to multiply by the orientation')
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MSH_RESTRICT_TO_PATCH: extracts the fields corresponding to the selected elements of a given patch, 
%       from the ones of a mesh struct, computed with msh_multipatch/msh_evaluate_element_list.
% The result is the same as calling msh_cartesian/msh_evaluate_element_list, but avoids recomputing.
%
function msh_ptc = msh_restrict_to_patch (msh, patch)

  msh_ptc.ndim = msh.ndim;
  msh_ptc.rdim = msh.rdim;
  msh_ptc.elem_list = msh.elem_list_of_patch{patch}(:).';

  msh_ptc.nel = msh.nel_per_patch(patch);
  msh_ptc.nel_dir = msh.nel_dir_of_patch{patch};
  msh_ptc.nqn = msh.nqn;
  msh_ptc.nqn_dir = msh.nqn_dir;

  if (isempty (msh_ptc.elem_list)), return, end
  
  Nelem = cumsum ([0, msh.nel_per_patch]);
  
  global_elem_list = Nelem(patch)+1:Nelem(patch+1);
  msh_ptc.quad_weights = msh.quad_weights(:,global_elem_list);
  msh_ptc.geo_map      = msh.geo_map(:,:,global_elem_list);
  msh_ptc.geo_map_jac  = msh.geo_map_jac(:,:,:,global_elem_list);
  msh_ptc.geo_map_der2 = msh.geo_map_der2(:,:,:,:,global_elem_list);
  msh_ptc.jacdet       = msh.jacdet(:,global_elem_list);
  msh_ptc.element_size = msh.element_size(:,global_elem_list);
  
end