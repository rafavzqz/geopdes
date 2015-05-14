% SP_PRECOMPUTE_PARAM: precompute all the fields, as in the space structure of the technical report, before mapping to the physical domain.
%  This function is used in vectorial spaces, before applying the map.
%
%     space = sp_precompute_param (space, msh);
%     space = sp_precompute_param (space, msh, 'option');
%
% INPUT:
%     
%    space: object representing the discrete function space (see sp_bspline).
%    msh: mesh object containing the quadrature information (see msh_cartesian)
%    'option', value: additional optional parameters, available options are:
%        nsh, connectivity, value (shape_functions), gradient (shape_function_gradients).
%     The value must be true or false. All the values are false by default.
%
% OUTPUT:
%
%    space: object containing the information of the input object, plus the 
%            fields of the old structure, that are listed below. If no option
%            is given all the fields are computed. If an option is given,
%            only the selected fields will be computed.
%
%    FIELD_NAME      (SIZE)                             DESCRIPTION
%    nsh             (1 x msh.nel vector)               actual number of shape functions per each element
%    connectivity    (nsh_max x msh.nel vector)         indices of basis functions that do not vanish in each element
%    shape_functions (msh.nqn x nsh_max x msh.nel)      basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%          (ndim x msh.nqn x nsh_max x msh.nel)         basis function gradients evaluated at each quadrature node in each element
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2015 Rafael Vazquez
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

function sp = sp_precompute_param (sp, msh, varargin)

  if (isempty (varargin))
    nsh = true;
    connectivity = true;
    value = true;
    gradient = true;
    hessian = true;
  else
    if (~rem (length (varargin), 2) == 0)
      error ('sp_precompute: options must be passed in the [option, value] format');
    end
    nsh = false;
    connectivity = false;
    value = false;
    gradient = false;
    hessian = false;
    for ii=1:2:length(varargin)-1
      if (strcmpi (varargin{ii}, 'connectivity'))
        connectivity = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'nsh'))
        nsh = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'value'))
        value = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'gradient'))
        gradient = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'hessian'))
        hessian = varargin{ii+1};
      else
        error ('sp_precompute_param: unknown option %s', varargin {ii});
      end
    end    
  end

  sp_univ = sp.sp_univ;

  for idim = 1:msh.ndim
    elem_list{idim} = 1:msh.nel_dir(idim);
  end
  
  if (nsh)
    for idim = 1:msh.ndim
      nsh_dim{idim} = sp_univ(idim).nsh(elem_list{idim});
    end

    [nsh_grid{1:msh.ndim}] = ndgrid (nsh_dim{:});
    nsh = 1;
    for idim = 1:msh.ndim
      nsh = nsh .* nsh_grid{idim};
    end
    sp.nsh = nsh(:)';
  end

  if (connectivity)
    for idim = 1:msh.ndim
      csize = ones (1, 2*msh.ndim);
      csize([idim, msh.ndim+idim]) = [sp_univ(idim).nsh_max, msh.nel_dir(idim)];
      crep = [sp_univ.nsh_max, msh.nel_dir];
      crep([idim, msh.ndim+idim]) = 1;

      conn{idim} = reshape (sp_univ(idim).connectivity(:,elem_list{idim}), csize);
      conn{idim} = repmat (conn{idim}, crep);
      conn{idim} = reshape (conn{idim}, [], msh.nel);
    end

    connectivity = zeros (sp.nsh_max, msh.nel);
    indices = ones (size (conn{1}));
    for idim = 1:msh.ndim
      indices = indices & conn{idim} ~= 0;
    end
    for idim = 1:msh.ndim
      conn{idim} = conn{idim}(indices);
    end
    connectivity(indices) = sub2ind ([sp.ndof_dir, 1], conn{:}); % The extra 1 makes things work in any dimension
    sp.connectivity = reshape (connectivity, sp.nsh_max, msh.nel);
    clear conn csize crep indices connectivity
  end

  if (value || gradient || hessian)
    shp = cell(1,msh.ndim); shg = cell(1,msh.ndim); shh = cell(1,msh.ndim);
    for idim = 1:msh.ndim
      ssize = ones (1, 3*msh.ndim);
      ssize([idim, msh.ndim+idim, 2*msh.ndim+idim]) = [msh.nqn_dir(idim), sp_univ(idim).nsh_max, msh.nel_dir(idim)];
      srep = [msh.nqn_dir, sp_univ.nsh_max, msh.nel_dir];
      srep([idim, msh.ndim+idim, 2*msh.ndim+idim]) = 1;
      shp{idim} = reshape (sp_univ(idim).shape_functions(:,:,elem_list{idim}), ssize);
      shp{idim} = repmat (shp{idim}, srep);
      shp{idim} = reshape (shp{idim}, msh.nqn, sp.nsh_max, msh.nel);
      shg{idim} = reshape (sp_univ(idim).shape_function_gradients(:,:,elem_list{idim}), ssize);
      shg{idim} = repmat (shg{idim}, srep);
      shg{idim} = reshape (shg{idim}, msh.nqn, sp.nsh_max, msh.nel);
      shh{idim} = reshape (sp_univ(idim).shape_function_hessians(:,:,elem_list{idim}), ssize);
      shh{idim} = repmat (shh{idim}, srep);
      shh{idim} = reshape (shh{idim}, msh.nqn, sp.nsh_max, msh.nel);
    end
    
    if (value)
      sp.shape_functions = 1;
      for idim = 1:msh.ndim
        sp.shape_functions = sp.shape_functions .* shp{idim};
      end
    end

    if (gradient)
      for idim = 1:msh.ndim
        shape_fun_grad = shg{idim};
        for jdim = setdiff (1:msh.ndim, idim)
          shape_fun_grad = shape_fun_grad .* shp{jdim};
        end
        sp.shape_function_gradients(idim,:,:,:) = shape_fun_grad;
      end
      sp.shape_function_gradients = reshape (sp.shape_function_gradients, ...
                                    msh.ndim, msh.nqn, sp.nsh_max, msh.nel);
    end

% To be changed with isprop, as soon as classdef is working
    if (hessian && isfield (struct (msh), 'geo_map_der2'))
      for idim = 1:msh.ndim
        shape_fun_hess = shh{idim};
        for jdim = setdiff (1:msh.ndim, idim)
          shape_fun_hess = shape_fun_hess .* shp{jdim};
        end
        sp.shape_function_hessians(idim,idim,:,:,:) = shape_fun_hess;
    
        for jdim = setdiff (1:msh.ndim, idim)
          shape_fun_hess = shg{idim} .* shg{jdim};
          for kdim = setdiff (1:msh.ndim, [idim, jdim])
            shape_fun_hess = shape_fun_hess .* shp{kdim};
          end
          sp.shape_function_hessians(idim,jdim,:,:,:) = shape_fun_hess;
        end
      end
    end
    
  end

end
