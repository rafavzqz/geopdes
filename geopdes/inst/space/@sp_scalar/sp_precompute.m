% SP_PRECOMPUTE: compute the basis functions in the whole mesh. The output
%  can consume a lot of memory. Use with care.
%
%     space = sp_precompute (space, msh);
%     space = sp_precompute (space, msh, 'option', value);
%
% INPUT:
%     
%    space: object representing the discrete function space (see sp_scalar).
%    msh: mesh object containing the quadrature information (see msh_cartesian)
%   'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
%            gradient   |      false      |  compute shape_function_gradients
%            hessian    |      false      |  compute shape_function_hessians
%            laplacian  |      false      |  compute shape_function_laplacians
%
% OUTPUT:
%
%     space: object containing the information of the input object, plus the 
%            fields of the old structure, that are listed below. If no option
%            is given all the fields are computed. If an option is given,
%            only the selected fields will be computed.
%
%    FIELD_NAME      (SIZE)                             DESCRIPTION
%    ncomp           (scalar)                           number of components of the functions of the space (actually, 1)
%    ndof            (scalar)                           total number of degrees of freedom
%    ndof_dir        (1 x ndim vector)                  degrees of freedom along each direction
%    nsh_max         (scalar)                           maximum number of shape functions per element
%    nsh             (1 x msh.nel vector)               actual number of shape functions per each element
%    connectivity    (nsh_max x msh.nel vector)         indices of basis functions that do not vanish in each element
%    shape_functions (msh.nqn x nsh_max x msh.nel)      basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%          (rdim x msh.nqn x nsh_max x msh.nel)         basis function gradients evaluated at each quadrature node in each element
%    shape_function_hessians
%          (rdim x rdim x msh.nqn x nsh_max x msh.nel)  basis function hessians evaluated at each quadrature node in each element
%    shape_function_laplacians 
%          (msh.nqn x nsh_max x msh.nel)                basis functions laplacians evaluated at each quadrature node in each element
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

  hessian_param = hessian || laplacian;
  grad_param = gradient || hessian_param;
  value_param = value || grad_param;
  
  sp_out = sp_precompute_param (sp, msh, 'value', value_param, 'gradient', grad_param, 'hessian', hessian_param);

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
