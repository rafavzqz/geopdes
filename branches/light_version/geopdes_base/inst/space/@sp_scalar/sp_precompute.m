% SP_PRECOMPUTE: precompute all the fields, as in the space structure of the technical report.
%
%     space = sp_precompute (space, msh);
%     space = sp_precompute (space, msh, 'option', value);
%
% INPUT:
%     
%    space: object representing the discrete function space (see sp_scalar).
%    msh: mesh object containing the quadrature information (see msh_cartesian)
%    'option', value: additional optional parameters, available options are:
%        nsh, connectivity, value (shape_functions), gradient (shape_function_gradients).
%     The value must be true or false. All the values are false by default.
%
% OUTPUT:
%
%     space: object containing the information of the input object, plus the 
%            fields of the old structure, that are listed below. If no option
%            is given all the fields are computed. If an option is given,
%            only the selected fields will be computed.
%
%    FIELD_NAME      (SIZE)                             DESCRIPTION
%    nsh             (1 x msh.nel vector)               actual number of shape functions per each element
%    connectivity    (nsh_max x msh.nel vector)         indices of basis functions that do not vanish in each element
%    shape_functions (msh.nqn x nsh_max x msh.nel)      basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%          (rdim x msh.nqn x nsh_max x msh.nel)         basis function gradients evaluated at each quadrature node in each element
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

function sp_out = sp_precompute (sp, msh, varargin)

  if (isempty (varargin))
    value = true;
    gradient = false;
    hessian = false;
    laplacian = false;
  else
    if (~rem (length (varargin), 2) == 0)
      error ('sp_precompute: options must be passed in the [option, value] format');
    end

    value = true;
    gradient = false;
    hessian = false;
    laplacian = false;
    for ii=1:2:length(varargin)-1
      if (strcmpi (varargin{ii}, 'value'))
        value = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'gradient'))
        gradient = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'hessian'))
        hessian = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'laplacian'))
        laplacian = varargin{ii+1};
      end
    end    
  end

  if (~isstruct (msh))
    msh = msh_precompute (msh);
  end
  
  if (isempty (varargin))
    sp_out = sp_precompute_param (sp, msh);
  else
    sp_out = sp_precompute_param (sp, msh, varargin{:});
  end

  switch (lower (sp.transform))
    case {'grad-preserving'}
      sp_out = sp_grad_preserving_transform (sp_out, msh, value, gradient, hessian, laplacian);
    case {'integral-preserving'}
      sp_out = sp_integral_preserving_transform (sp_out, msh, value);
      if (gradient || hessian || laplacian)
        warning ('Derivatives not implemented for integral-preserving transformation')
      end
  end

end
