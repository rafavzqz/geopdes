% SP_EVALUATE_ELEMENT_LIST_PARAM: compute the basis functions, in the parametric domain, in a given list of elements.
%
%     sp = sp_evaluate_element_list_param (space, msh_elems, 'option1', value1, ...)
%
% INPUTS:
%     
%    space:   object defining the space of discrete functions (see sp_vector)
%    msh_elems: msh structure containing the information of quadrature or
%               visualization points, for a given list of elements 
%               (see msh_cartesian/msh_evaluate_element_list)
%   'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
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
%    ndof_dir        (ncomp_param x ndim matrix)                for each component, number of degrees of freedom along each direction
%    nsh_max         (scalar)                                   maximum number of shape functions per element
%    nsh             (1 x msh_elems.nel vector)                 actual number of shape functions per each element
%    connectivity    (nsh_max x msh_elems.nel vector)           indices of basis functions that do not vanish in each element
%    shape_functions (ncomp_param x msh_elems.nqn x nsh_max x msh_elems.nel)  basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%       (ncomp_param x ndim x msh_elems.nqn x nsh_max x msh_elems.nel) basis function gradients evaluated at each quadrature node in each element
%    shape_function_hessians
%       (ncomp_param x ndim x ndim x msh_elems.nqn x nsh_max x msh_elems.nel) basis function hessians evaluated at each quadrature node in each element
%    shape_function_divs (msh_elems.nqn x nsh_max x msh_elems.nel)     basis function divergence evaluated at each quadrature node in each element
%    shape_function_curls 
%         2D:  (msh_elems.nqn x nsh_max x msh_elems.nel)               basis function curl evaluated at each quadrature node in each element
%         3D:  (3 x msh_elems.nqn x nsh_max x msh_elems.nel)        
%    shape_function_third_derivatives
%       (ncomp_param x ndim x ndim x ndim x msh_elems.nqn x nsh_max x msh_elems.nel) basis function third derivatives evaluated at each quadrature node in each element
%    shape_function_fourth_derivatives
%       (ncomp_param x ndim x ndim x ndim x rdim x msh_elems.nqn x nsh_max x msh_elems.nel) basis function fourth derivatives evaluated at each quadrature node in each element
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

function sp = sp_evaluate_element_list_param (space, msh, varargin)

value = true;
gradient = false;
divergence = false;
curl = false;
hessian = false;
third_derivative = false;
fourth_derivative = false;

if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_evaluate_element_list_param: options must be passed in the [option, value] format');
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

first_der = gradient || divergence || curl;
for icomp = 1:space.ncomp_param
  sp_col_scalar(icomp) = sp_evaluate_element_list_param (space.scalar_spaces{icomp}, msh, 'value', value, 'gradient', first_der, ...,
                                                'hessian', hessian, 'third_derivative', third_derivative, ...,
                                                'fourth_derivative', fourth_derivative);
end

ndof_scalar = [sp_col_scalar.ndof];

ndof = sum (ndof_scalar);
nsh  = zeros (1, msh.nel);
connectivity = [];
aux = 0;
for icomp = 1:space.ncomp_param
  ndof_dir(icomp,:) = sp_col_scalar(icomp).ndof_dir;
  nsh = nsh + sp_col_scalar(icomp).nsh(:)';
  
  inds = find (sp_col_scalar(icomp).connectivity);
  sp_col_scalar(icomp).connectivity(inds) = sp_col_scalar(icomp).connectivity(inds) + aux;
  connectivity = [connectivity; sp_col_scalar(icomp).connectivity];
  aux = aux + ndof_scalar(icomp);
end

sp = struct('nsh_max', space.nsh_max, 'nsh', nsh, 'ndof', ndof,  ...
            'ndof_dir', ndof_dir, 'connectivity', connectivity, ...
            'ncomp', space.ncomp, 'ncomp_param', space.ncomp_param);

if (value)
  sp.shape_functions = zeros (sp.ncomp_param, msh.nqn, sp.nsh_max, msh.nel);
  for icomp = 1:space.ncomp_param
    indices = space.cumsum_nsh(icomp)+(1:sp_col_scalar(icomp).nsh_max);
    sp.shape_functions(icomp,:,indices,:) = sp_col_scalar(icomp).shape_functions;
  end
end

if (gradient || curl || divergence)
  shape_fun_grads = zeros (space.ncomp_param, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);

  for icomp = 1:space.ncomp_param
    indices = space.cumsum_nsh(icomp)+(1:sp_col_scalar(icomp).nsh_max);
    shape_fun_grads(icomp,:,:,indices,:) = sp_col_scalar(icomp).shape_function_gradients;
  end
  
  if (gradient)
    sp.shape_function_gradients = shape_fun_grads;
  end

  if (divergence)
    sp.shape_function_divs = zeros (msh.nqn, sp.nsh_max, msh.nel);
    for icomp = 1:space.ncomp_param
      sp.shape_function_divs = sp.shape_function_divs + ...
        reshape (shape_fun_grads(icomp,icomp,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
    end
  end

  if (curl)
    if (space.ncomp_param == 2)
      sp.shape_function_curls = reshape (shape_fun_grads(2,1,:,:,:) - ...
	 	       shape_fun_grads(1,2,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
    elseif (space.ncomp_param == 3)
      shape_fun_curls = zeros (space.ncomp, msh.nqn, sp.nsh_max, msh.nel);
      for icomp = 1:space.ncomp
        ind1 = mod(icomp,3) + 1;
        ind2 = mod(ind1, 3) + 1;
        shape_fun_curls(icomp,:,:,:) = reshape (shape_fun_grads(ind2,ind1,:,:,:) - ...
                   shape_fun_grads(ind1,ind2,:,:,:), 1, msh.nqn, sp.nsh_max, msh.nel);
      end
      sp.shape_function_curls = shape_fun_curls;
    end
  end
end

if (hessian)
  sp.shape_function_hessians = zeros (space.ncomp, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);
  for icomp = 1:space.ncomp_param
    indices = space.cumsum_nsh(icomp)+(1:sp_col_scalar(icomp).nsh_max);
    sp.shape_function_hessians(icomp,:,:,:,indices,:) = sp_col_scalar(icomp).shape_function_hessians;
  end
end

if (third_derivative)
  sp.shape_function_third_derivatives = zeros (space.ncomp, msh.ndim, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);
  for icomp = 1:space.ncomp_param
    indices = space.cumsum_nsh(icomp)+(1:sp_col_scalar(icomp).nsh_max);
    sp.shape_function_third_derivatives(icomp,:,:,:,:,indices,:) = sp_col_scalar(icomp).shape_function_third_derivatives;
  end
end

if (fourth_derivative)
  sp.shape_function_fourth_derivatives = zeros (space.ncomp, msh.ndim, msh.ndim, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);
  for icomp = 1:space.ncomp_param
    indices = space.cumsum_nsh(icomp)+(1:sp_col_scalar(icomp).nsh_max);
    sp.shape_function_fourth_derivatives(icomp,:,:,:,:,:,indices,:) = sp_col_scalar(icomp).shape_function_fourth_derivatives;
  end
end

end
