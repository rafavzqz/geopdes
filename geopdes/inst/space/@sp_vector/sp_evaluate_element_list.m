% SP_EVALUATE_ELEMENT_LIST: compute the basis functions in a given list of elements.
%
%     sp = sp_evaluate_element_list (space, msh_elems, 'option1', value1, ...)
%
% INPUTS:
%     
%    space:     object defining the space of discrete functions (see sp_vector)
%    msh_elems: msh structure containing the information of quadrature or
%               visualization points, for a given list of elements 
%               (see msh_cartesian/msh_evaluate_element_list)
%   'option', value: additional optional parameters, currently available options are:
%            
%              Name             |   Default value |  Meaning
%           --------------------+-----------------+----------------------------------
%            value              |      true       |  compute shape_functions
%            gradient           |      false      |  compute shape_function_gradients
%            divergence         |      false      |  compute shape_function_divs
%            curl               |      false      |  compute shape_function_curls
%            third_derivative   |      false      |  compute shape_function_third_derivatives
%            fourth_derivative  |      false      |  compute shape_function_fourth_derivatives
%
% OUTPUT:
%
%    sp: struct representing the discrete function space, with the following fields:
%              (see the article for a detailed description)
%
%    FIELD_NAME      (SIZE)                                     DESCRIPTION
%    ncomp           (scalar)                                   number of components of the functions of the space
%    ndof            (scalar)                                   total number of degrees of freedom
%    ndof_dir        (ncomp x ndim matrix)                      for each component, number of degrees of freedom along each direction
%    nsh_max         (scalar)                                   maximum number of shape functions per element
%    nsh             (1 x msh_elems.nel vector)                 actual number of shape functions per each element
%    connectivity    (nsh_max x msh_elems.nel vector)           indices of basis functions that do not vanish in each element
%    shape_functions (ncomp x msh_elems.nqn x nsh_max x msh_elems.nel)  basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%       (ncomp x rdim x msh_elems.nqn x nsh_max x msh_elems.nel)    basis function gradients evaluated at each quadrature node in each element
%    shape_function_divs (msh_elems.nqn x nsh_max x msh_elems.nel)  basis function divergence evaluated at each quadrature node in each element
%    shape_function_curls 
%         2D:  (msh_elems.nqn x nsh_max x msh_elems.nel)            basis function curl evaluated at each quadrature node in each element
%         3D:  (3 x msh_elems.nqn x nsh_max x msh_elems.nel)        
%    shape_function_hessians
%        rdim x rdim x msh_elems.nqn x nsh_max x msh_elems.nel) basis function hessians evaluated at each quadrature node in each element
%    shape_function_third_derivatives
%       (rdim x rdim x rdim x msh_elems.nqn x nsh_max x msh_elems.nel) basis function third derivatives evaluated at each quadrature node in each element
%    shape_function_fourth_derivatives
%       (rdim x rdim x rdim x rdim x msh_elems.nqn x nsh_max x msh_elems.nel) basis function fourth derivatives evaluated at each quadrature node in each element
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011, 2015, 2019 Rafael Vazquez
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

function sp = sp_evaluate_element_list (space, msh, varargin)

value = true;
gradient = false;
divergence = false;
curl = false;
hessian = false;
third_derivative = false;
fourth_derivative = false;

if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_evaluate_element_list: options must be passed in the [option, value] format');
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

fourth_param = fourth_derivative;
third_param = third_derivative || fourth_param;
hessian_param = hessian || third_param;
grad_param = gradient || hessian_param || divergence || curl;
value_param = value || grad_param;
div_param = false; curl_param = false;
switch (lower (space.transform))
  case {'curl-preserving'}
    curl_param = curl;
  case {'div-preserving'}
    div_param = divergence;
end

sp = sp_evaluate_element_list_param (space, msh, 'value', value_param, 'gradient', grad_param, 'divergence', div_param,...
                                     'curl', curl_param, 'hessian', hessian, 'third_derivative', third_param, ...
                                     'fourth_derivative', fourth_param);

if (isempty (msh.elem_list))
  return
end

switch (lower (space.transform))
  case {'grad-preserving'}
    sp = sp_vector_grad_preserving_transform (sp, msh, value, gradient, curl, divergence, hessian, third_derivative, fourth_derivative);
  case {'curl-preserving'}
    sp = sp_vector_curl_preserving_transform (sp, msh, value, curl);
    if (gradient || divergence || hessian || third_derivative || fourth_derivative)
      warning ('Gradient, divergence, hessian, and third and fourth derivatives not implemented for curl-preserving transformation')
    end
  case {'div-preserving'}
    sp = sp_vector_div_preserving_transform (sp, msh, value, gradient, curl, divergence);
    if (hessian || third_derivative || fourth_derivative)
      warning ('Hessian and third and fourth derivatives not implemented for div-preserving transformation')
    end
end

end
