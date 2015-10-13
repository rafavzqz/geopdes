% SP_PRECOMPUTE: precompute all the fields, as in the space structure of the technical report.
%
%     space = sp_precompute (space, msh);
%     space = sp_precompute (space, msh, 'option', value);
%
% INPUT:
%     
%    space: object representing the discrete function space (see sp_vector_curl_transform).
%    msh: mesh object containing the quadrature information (see msh_cartesian)
%    'option', value: additional optional parameters, available options are:
%        nsh, connectivity, value (shape_functions), curl (shape_function_curls)
%     The value must be true or false. All the values are false by default.
%
% OUTPUT:
%
%    space: object containing the information of the input object, plus the 
%            fields of the old structure, that are listed below. If no option
%            is given all the fields are computed. If an option is given,
%            only the selected fields will be computed.
%
%    FIELD_NAME      (SIZE)                                    DESCRIPTION
%    nsh             (1 x msh.nel vector)                      actual number of shape functions per each element
%    connectivity    (nsh_max x msh.nel vector)                indices of basis functions that do not vanish in each element
%    shape_functions (rdim x msh.nqn x nsh_max x msh.nel)      basis functions evaluated at each quadrature node in each element
%  In 2D
%    shape_function_curls (msh.nqn x nsh_max x msh.nel)        basis function curls evaluated at each quadrature node in each element
%  In 3D
%    shape_function_curls (rdim x msh.nqn x nsh_max x msh.nel) basis function curls evaluated at each quadrature node in each element
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2010, 2011, 2015 Rafael Vazquez
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

function sp = sp_precompute (sp, msh, varargin)

  if (isempty (varargin))
    nsh = true;
    connectivity = true;
    value = true;
    curl = false;
    scalar_gradient = false;
  else
    if (~rem (length (varargin), 2) == 0)
      error ('sp_precompute: options must be passed in the [option, value] format');
    end

    nsh = false;
    connectivity = false;
    value = false;
    curl = false;
    scalar_gradient = false;
    for ii=1:2:length(varargin)-1
      if (strcmpi (varargin{ii}, 'connectivity'))
        connectivity = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'nsh'))
        nsh = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'value'))
        value = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'curl'))
        curl = varargin{ii+1};
      else
        error ('sp_precompute: unknown option %s', varargin {ii});
      end
    end
    if (curl)
      scalar_gradient = true;
    end
  end

% Precompute everything for each component in the parametric domain
% These structures will not be stored in memory
  scalar_sp = cell (sp.ncomp_param,1);
  for icomp = 1:sp.ncomp_param
    scalar_sp{icomp} = sp_precompute (sp.scalar_spaces{icomp}, msh, 'nsh', nsh, 'connectivity', ...
                    connectivity, 'value', value, 'gradient', scalar_gradient);
  end
  
  if (nsh)
    sp.nsh  = zeros (1, msh.nel);
    for icomp = 1:sp.ncomp_param
      sp.nsh = sp.nsh + scalar_sp{icomp}.nsh(:)';
    end
  end

  if (connectivity)
    sp.connectivity = [];
    for icomp = 1:sp.ncomp_param
      sp.connectivity = [sp.connectivity; scalar_sp{icomp}.connectivity+sp.cumsum_ndof(icomp)];
    end
  end

% From here we apply the transformation
  if (value || curl)
    if (isempty (msh.geo_map_jac))
      msh = msh_precompute (msh, 'geo_map_jac', true);
    end
    [JinvT, jacdet] = geopdes_invT__ (msh.geo_map_jac);
  end

  if (value)
    JinvT = reshape (JinvT, [msh.rdim, msh.ndim, msh.nqn, msh.nel]);
    shape_functions = zeros (msh.ndim, msh.nqn, sp.nsh_max, msh.nel);
    for icomp = 1:sp.ncomp_param
      indices = sp.cumsum_nsh(icomp)+(1:scalar_sp{icomp}.nsh_max);
      shape_functions(icomp,:,indices,:) = scalar_sp{icomp}.shape_functions;
    end
    sp.shape_functions = geopdes_prod__ (JinvT, shape_functions);
    clear shape_functions
  end
  
  if (curl)
    if (sp.ncomp_param == 2)
      shape_function_curls = zeros (msh.nqn, sp.nsh_max, msh.nel);
      for icomp = 1:sp.ncomp_param
        ind = setdiff (1:msh.ndim, icomp);
        indices = sp.cumsum_nsh(icomp)+(1:scalar_sp{icomp}.nsh_max);

        shape_function_curls(:,indices,:) = (-1)^icomp * ...
          reshape (scalar_sp{icomp}.shape_function_gradients(ind,:,:,:), msh.nqn, scalar_sp{icomp}.nsh_max, msh.nel);
      end
      jacdet = reshape (jacdet, msh.nqn, 1, msh.nel);
      sp.shape_function_curls = bsxfun (@rdivide, shape_function_curls, jacdet);
  
    elseif (sp.ncomp_param == 3)
      shape_function_curls = zeros (sp.ncomp_param, msh.nqn, sp.nsh_max, msh.nel);
      for icomp = 1:sp.ncomp_param
        ind(1) = mod (icomp, 3) + 1; ind(2) = mod (ind(1), 3) + 1;
        ind2 = fliplr (ind);
        indices = sp.cumsum_nsh(icomp)+(1:scalar_sp{icomp}.nsh_max);
        for ii = 1:numel(ind)
          shape_function_curls(ind(ii),:,indices,:) = (-1)^(ii-1) * ...
            reshape (scalar_sp{icomp}.shape_function_gradients(ind2(ii),:,:,:), msh.nqn, scalar_sp{icomp}.nsh_max, msh.nel);
        end
      end
    
      jacdet = reshape (jacdet, 1, msh.nqn, 1, msh.nel);
      shape_function_curls = geopdes_prod__ (msh.geo_map_jac, shape_function_curls);
      sp.shape_function_curls = bsxfun (@rdivide, shape_function_curls, jacdet);
    end
  end

end
