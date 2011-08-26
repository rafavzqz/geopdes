% SP_PRECOMPUTE: precompute all the fields, as in the space structure of the technical report.
%
%     space = sp_precompute (space, msh);
%     space = sp_precompute (space, msh, 'option', value);
%
% INPUT:
%     
%    space: object representing the discrete function space (see sp_vector_2d_curl_transform).
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
%    FIELD_NAME      (SIZE)                                  DESCRIPTION
%    nsh             (1 x msh.nel vector)                    actual number of shape functions per each element
%    connectivity    (nsh_max x msh.nel vector)              indices of basis functions that do not vanish in each element
%    shape_functions (2 x msh.nqn x nsh_max x msh_col.nel)   basis functions evaluated at each quadrature node in each element
%    shape_function_curls (msh.nqn x nsh_max x msh.nel)      basis function curl evaluated at each quadrature node in each element
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2010, 2011 Rafael Vazquez
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
    curl = true;
    scalar_gradient = true;
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
  sp1 = sp_precompute_param (sp.sp1, msh, 'nsh', nsh, 'connectivity', ...
                    connectivity, 'value', value, 'gradient', scalar_gradient);
  sp2 = sp_precompute_param (sp.sp2, msh, 'nsh', nsh, 'connectivity', ...
                    connectivity, 'value', value, 'gradient', scalar_gradient);

  if (nsh)
    sp.nsh = sp1.nsh(:)' + sp2.nsh(:)';
  end

  if (connectivity)
    sp.connectivity = [sp1.connectivity; sp2.connectivity+sp1.ndof];
  end

% From here we apply the transformation
  if (value || curl)
    if (isempty (msh.geo_map_jac))
      msh = msh_precompute (msh, 'geo_map_jac', true);
    end
    [JinvT, jacdet] = geopdes_invT__ (msh.geo_map_jac);
  end

  if (value)
    JinvT = reshape (JinvT, [2, 2, msh.nqn, msh.nel]);
    sp.shape_functions = zeros (2, msh.nqn, sp.nsh_max, msh.nel);
    sp.shape_functions(1,:,1:sp1.nsh_max,:)            = sp1.shape_functions;
    sp.shape_functions(2,:,sp1.nsh_max+1:sp.nsh_max,:) = sp2.shape_functions;
    sp.shape_functions = geopdes_prod__ (JinvT, sp.shape_functions);
  end

  if (curl)
    sp.shape_function_curls = zeros (msh.nqn, sp.nsh_max, msh.nel);
    sp.shape_function_curls(:, 1:sp1.nsh_max, :) = reshape ...
      (-sp1.shape_function_gradients(2,:,:,:), msh.nqn, sp1.nsh_max, msh.nel);
    sp.shape_function_curls(:, sp1.nsh_max+1:sp.nsh_max, :) = reshape ...
       (sp2.shape_function_gradients(1,:,:,:), msh.nqn, sp2.nsh_max, msh.nel);
    for ii=1:sp.nsh_max
      sp.shape_function_curls(:,ii,:) = ...
         reshape (sp.shape_function_curls(:,ii,:), msh.nqn, msh.nel)./jacdet;
    end
  end


end
