% SP_PRECOMPUTE: compute the basis functions in the whole mesh. The output
%  can consume a lot of memory. Use with care.
%
%     space = sp_precompute (space, msh);
%     space = sp_precompute (space, msh, 'option', value);
%
% INPUT:
%     
%    space: object representing the discrete function space (see sp_vector)
%    msh: mesh object containing the quadrature information (see msh_cartesian)
%   'option', value: additional optional parameters, currently available options are:
%            
%              Name             |   Default value |  Meaning
%           --------------------+-----------------+----------------------------------
%            value              |      true       |  compute shape_functions
%            gradient           |      false      |  compute shape_function_gradients
%            divergence         |      false      |  compute shape_function_divs
%            curl               |      false      |  compute shape_function_curls
%            hessian            |      false      |  compute shape_function_hessians
%            third_derivative   |      false      |  compute shape_function_third_derivatives
%            fourth_derivative  |      false      |  compute shape_function_fourth_derivatives
%
% OUTPUT:
%
%    space: object containing the information of the input object, plus the 
%            fields of the old structure, that are listed below. If no option
%            is given all the fields are computed. If an option is given,
%            only the selected fields will be computed.
%
%    FIELD_NAME      (SIZE)                               DESCRIPTION
%    ncomp           (scalar)                              number of components of the functions of the space
%    ndof            (scalar)                              total number of degrees of freedom
%    ndof_dir        (ncomp x ndim matrix)                 for each component, number of degrees of freedom along each direction
%    nsh_max         (scalar)                              maximum number of shape functions per element
%    nsh             (1 x msh.nel vector)                  actual number of shape functions per each element
%    connectivity    (nsh_max x msh.nel vector)            indices of basis functions that do not vanish in each element
%    shape_functions (ncomp x msh.nqn x nsh_max x msh.nel) basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%             (ncomp x rdim x msh.nqn x nsh_max x msh.nel) basis function gradients evaluated at each quadrature node in each element
%    shape_function_divs (msh.nqn x nsh_max x msh.nel)     basis function divergence evaluated at each quadrature node in each element
%    shape_function_curls
%         2D:  (msh.nqn x nsh_max x msh.nel)               basis function curl evaluated at each quadrature node in each element
%         3D:  (3 x msh.nqn x nsh_max x msh.nel)        
%    shape_function_hessians
%          (ncomp x rdim x rdim x msh.nqn x nsh_max x msh.nel)  basis function hessians evaluated at each quadrature node in each element
%    shape_function_third_derivatives
%       (ncomp x rdim x rdim x rdim x msh.nqn x nsh_max x msh.nel) basis function third derivatives evaluated at each quadrature node in each element
%    shape_function_fourth_derivatives
%       (ncomp x rdim x rdim x rdim x rdim x msh.nqn x nsh_max x msh.nel) basis function fourth derivatives evaluated at each quadrature node in each element

%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2015, 2019 Rafael Vazquez
% Copyright (C) 2023 Pablo Antolin, Luca Coradello
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

value = true;
gradient = false;
divergence = false;
curl = false;
hessian = false;
third_derivative = false;
fourth_derivative = false;

if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_precompute: options must be passed in the [option, value] format');
  end
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'value'))
      value = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'gradient'))
      gradient = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'curl'))
      curl = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'divergence'))
      divergence = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'hessian'))
      hessian = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'third_derivative'))
      third_derivative = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'fourth_derivative'))
      fourth_derivative = varargin {ii+1};
    else
      warning ('Ignoring unknown option %s', varargin {ii});
    end
  end
end

if (~isstruct (msh))
  msh = msh_precompute (msh);
end


fourth_param = fourth_derivative;
third_param = third_derivative || fourth_param;
hessian_param = hessian || third_param;
grad_param = gradient || hessian_param || divergence || curl;
value_param = value || grad_param;
div_param = false; curl_param = false;

switch (lower (sp.transform))
  case {'curl-preserving'}
    curl_param = curl;
  case {'div-preserving'}
    div_param = divergence;
end

sp_out = sp_precompute_param (sp, msh, 'value', value_param, 'gradient', grad_param, 'divergence', div_param, ...
                              'curl', curl_param, 'hessian', hessian, 'third_derivative', third_param, ...
                              'fourth_derivative', fourth_param);


switch (lower (sp.transform))
  case {'grad-preserving'}
    sp_out = sp_vector_grad_preserving_transform (sp_out, msh, value, gradient, curl, divergence, hessian, third_derivative, fourth_derivative);
  case {'curl-preserving'}
    sp_out = sp_vector_curl_preserving_transform (sp_out, msh, value, curl);
    if (gradient || divergence || hessian || third_derivative || fourth_derivative)
      warning ('Gradient, divergence, hessian, and third and fourth derivatives not implemented for curl-preserving transformation')
    end
  case {'div-preserving'}
    sp_out = sp_vector_div_preserving_transform (sp_out, msh, value, gradient, curl, divergence);
    if (hessian || third_derivative || fourth_derivative)
      warning ('Hessian and third and fourth derivatives not implemented for div-preserving transformation')
    end
end

end
