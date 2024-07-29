% SP_EVALUATE_ELEMENT_LIST: compute the basis functions in a given list of elements.
%
%     sp = sp_evaluate_element_list (space, msh_elems, 'option1', value1, ...)
%
% INPUTS:
%     
%    space:     object defining the space of discrete functions (see sp_multipatch_C1)
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
% Copyright (C) 2015, 2017, 2022 Rafael Vazquez
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

  is_scalar = isa (space.sp_patch{1}, 'sp_scalar');

  sp.npatch = space.npatch;
  sp.ncomp = space.ncomp;
  sp.ndof = space.ndof;

  if (isempty (msh.elem_list)), return, end

  sp.nsh_max = [];
  if (is_scalar)
    fields = {'nsh', 'connectivity', 'shape_functions', 'shape_function_gradients', 'shape_function_hessians', 'shape_function_laplacians'};
    cat_position = [2, 2, 3, 4, 5, 3];
%   else
%     fields = {'nsh', 'connectivity', 'shape_functions', 'shape_function_gradients', 'shape_function_hessians', ...
%       'shape_function_laplacians', 'shape_function_divs', 'shape_function_curls'};
%     if (msh.ndim == 2 || msh.ndim == 1)
%       cat_position = [2, 2, 4, 5, 6, 4, 3, 3];
%     elseif (msh.ndim == 3)
%       cat_position = [2, 2, 4, 5, 6, 4, 3, 4];
%     end
  end
  for ii = 1:numel(fields)
    sp.(fields{ii}) = [];
  end
  
  for iptc = 1:space.npatch
%     msh_patch = msh_evaluate_element_list (msh.msh_patch{iptc}, msh.elem_list_of_patch{iptc}); % Not working anymore
    msh_patch = msh_restrict_to_patch (msh, iptc);

    sp_patch = sp_evaluate_element_list (space.sp_patch{iptc}, msh_patch, varargin{:});
    [Cpatch, Cpatch_cols] = sp_compute_Cpatch (space, iptc);
    
    nsh = zeros (1, msh_patch.nel);
    connectivity = zeros (sp_patch.nsh_max, msh_patch.nel);
    shape_funs = zeros (msh_patch.nqn, sp_patch.nsh_max, msh_patch.nel);
    shape_fun_grads = zeros (msh.rdim, msh_patch.nqn, sp_patch.nsh_max, msh_patch.nel);
    shape_fun_hess = zeros (msh.rdim, msh.rdim, msh_patch.nqn, sp_patch.nsh_max, msh_patch.nel);
    shape_fun_lapl = zeros (msh_patch.nqn, sp_patch.nsh_max, msh_patch.nel);
    
    for iel = 1:msh_patch.nel
      conn_iel = sp_patch.connectivity(:,iel);
      [~,jj] = find (Cpatch(conn_iel,:));
      col_indices = unique(jj);
      funs = Cpatch_cols(col_indices);
      nsh(iel) = numel (funs);
      connectivity(1:nsh(iel),iel) = funs;
      Cpatch_iel = Cpatch(conn_iel, col_indices);

      if (isfield (sp_patch, 'shape_functions'))
        shape_funs(:,1:nsh(iel),iel) = sp_patch.shape_functions(:,:,iel) * Cpatch_iel;
      end
      if (isfield (sp_patch, 'shape_function_gradients'))
        for idim = 1:msh.rdim
          shape_fun_grads(idim,:,1:nsh(iel),iel) = ...
            reshape (sp_patch.shape_function_gradients(idim,:,:,iel), msh_patch.nqn, sp_patch.nsh_max) * Cpatch_iel;
        end
      end
      if (isfield (sp_patch, 'shape_function_laplacians'))
        shape_fun_lapl(:,1:nsh(iel),iel) = sp_patch.shape_function_laplacians(:,:,iel) * Cpatch_iel;
      end
      if (isfield (sp_patch, 'shape_function_hessians'))
        for idim = 1:msh.rdim
          for jdim = 1:msh.rdim
          shape_fun_hess(idim,jdim,:,1:nsh(iel),iel) = ...
            reshape (sp_patch.shape_function_hessians(idim,jdim,:,:,iel), msh_patch.nqn, sp_patch.nsh_max) * Cpatch_iel;
          end
        end
      end
    end

    nsh_max = max(sp.nsh); nsh_max_patch = max(nsh);
    sp.nsh = cat (2, sp.nsh, nsh);
    if (msh_patch.nel > 0)
      if (nsh_max > nsh_max_patch)
        connectivity(nsh_max_patch+1:nsh_max,:) = 0;
        if (isfield (sp_patch, 'shape_functions')); shape_funs(:,nsh_max_patch+1:nsh_max,:) = 0; end
        if (isfield (sp_patch, 'shape_function_gradients')); shape_fun_grads(:,:,nsh_max_patch+1:nsh_max,:) = 0; end
        if (isfield (sp_patch, 'shape_function_hessians')); shape_fun_hess(:,:,:,nsh_max_patch+1:nsh_max,:) = 0; end
        if (isfield (sp_patch, 'shape_function_laplacians')); shape_fun_lapl(:,nsh_max_patch+1:nsh_max,:) = 0; end
      elseif (nsh_max_patch > nsh_max)
        sp.connectivity(nsh_max+1:nsh_max_patch,:) = 0;
        if (isfield (sp_patch, 'shape_functions')); sp.shape_functions(:,nsh_max+1:nsh_max_patch,:) = 0; end
        if (isfield (sp_patch, 'shape_function_gradients')); sp.shape_function_gradients(:,:,nsh_max+1:nsh_max_patch,:) = 0; end
        if (isfield (sp_patch, 'shape_function_hessians')); sp.shape_function_hessians(:,:,:,nsh_max+1:nsh_max_patch,:) = 0; end
        if (isfield (sp_patch, 'shape_function_laplacians')); sp.shape_function_laplacians(:,nsh_max+1:nsh_max_patch,:) = 0; end
      end
      sp.connectivity = cat (2, sp.connectivity, connectivity);
      if (isfield (sp_patch, 'shape_functions')); sp.shape_functions = cat (3, sp.shape_functions, shape_funs); end
      if (isfield (sp_patch, 'shape_function_gradients')); sp.shape_function_gradients = cat (4, sp.shape_function_gradients, shape_fun_grads); end
      if (isfield (sp_patch, 'shape_function_hessians')); sp.shape_function_hessians = cat (5, sp.shape_function_hessians, shape_fun_hess); end
      if (isfield (sp_patch, 'shape_function_laplacians')); sp.shape_function_laplacians = cat (3, sp.shape_function_laplacians, shape_fun_lapl); end
    end
    
%     for ii = 1:numel(fields)
%       if (isfield (sp_patch, fields{ii}))
%         sp.(fields{ii}) = cat (cat_position(ii), sp.(fields{ii}), sp_patch.(fields{ii}));
%       end
%     end
  end
  sp.nsh_max = max (sp.nsh);

  for ii = 1:numel(fields)
    if (isempty (sp.(fields{ii})))
      sp = rmfield (sp, fields{ii});
    elseif (ii > 1)
      field_size = size (sp.(fields{ii}));
      inds = repmat ({':'}, 1, numel(field_size));
      inds{cat_position(ii)-1} = 1:sp.nsh_max;
      sp.(fields{ii}) = sp.(fields{ii})(inds{:});
    end
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