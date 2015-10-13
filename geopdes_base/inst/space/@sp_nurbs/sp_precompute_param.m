% SP_PRECOMPUTE_PARAM: precompute all the fields, as in the space structure of the technical report, before mapping to the physical domain.
%  This function is used in vectorial spaces, before applying the map.
%
%     space = sp_precompute_param (space, msh); computes all the fields
%     space = sp_precompute_param (space, msh, 'option');  only computes the selected fields
%
% INPUT:
%     
%    space: object representing the discrete function space (see sp_nurbs).
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

function sp = sp_precompute_param (space, msh, varargin)

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

% For NURBS, low order derivatives are used to compute high order derivatives
% They are removed at the end of the function
% The connectivity is also needed in bsp_2_nrb__
  conn_param = connectivity;
  value_param = value;
  gradient_param = gradient;
  if (hessian)
    gradient_param = true;
    value_param = true;
    conn_param = true;
  elseif (gradient)
    value_param = true;
    conn_param = true;
  elseif (value)
    conn_param = true;
  end

%   sp = sp_precompute_param (space.spline_space, msh, varargin{:});
  sp = sp_precompute_param (space.spline_space, msh, 'nsh', nsh, 'connectivity', conn_param, ...
         'value', value_param, 'gradient', gradient_param, 'hessian', hessian);

  if (value || gradient || hessian)
    sp = bsp_2_nrb__ (sp, msh, space.weights);
    
    if (~value)
      sp.shape_functions = [];
    end
    if (~gradient)
      sp.shape_function_gradients = [];
    end
    if (~connectivity)
      sp.connectivity = [];
    end
  end

  
  
end
