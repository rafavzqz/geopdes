% SP_EVALUATE_ELEMENT_LIST: compute the basis functions in a given list of elements.
%
%     sp = sp_evaluate_element_list (space, msh_elems, 'option1', value1, ...)
%
% INPUTS:
%     
%    space:     object defining the space of discrete functions (see sp_bspline)
%    msh_elems: msh structure containing the information of quadrature or
%               visualization points, for a given list of elements 
%               (see msh_cartesian/msh_evaluate_element_list)
%   'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
%            gradient   |      false      |  compute shape_function_gradients
%            hessian    |      false      |  compute shape_function_hessians
%
% OUTPUT:
%
%    sp: struct representing the discrete function space, with the following fields:
%              (see the article for a detailed description)
%
%    FIELD_NAME      (SIZE)                                     DESCRIPTION
%    ncomp           (scalar)                                   number of components of the functions of the space (actually, 1)
%    ndof            (scalar)                                   total number of degrees of freedom
%    ndof_dir        (1 x ndim vector)                          degrees of freedom along each direction
%    nsh_max         (scalar)                                   maximum number of shape functions per element
%    nsh             (1 x msh_elems.nel vector)                 actual number of shape functions per each element
%    connectivity    (nsh_max x msh_elems.nel vector)           indices of basis functions that do not vanish in each element
%    shape_functions (msh_elems.nqn x nsh_max x msh_elems.nel)  basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%       (rdim x msh_elems.nqn x nsh_max x msh_elems.nel)        basis function gradients evaluated at each quadrature node in each element
%    shape_function_hessians
%        rdim x rdim x msh_elems.nqn x nsh_max x msh_elems.nel) basis function hessians evaluated at each quadrature node in each element
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
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

function sp = sp_evaluate_element_list (space, msh, varargin)

value = true;
gradient = false;
grad_param = false;
hessian = false;
if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_evaluate_element_list: options must be passed in the [option, value] format');
  end
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'value'))
      value = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'gradient'))
      gradient = varargin {ii+1};
      grad_param = gradient;
    elseif (strcmpi (varargin {ii}, 'hessian'))
      hessian = varargin {ii+1};
    else
      error ('sp_evaluate_element_list: unknown option %s', varargin {ii});
    end
  end
end

if (hessian)
  grad_param = true;
end
sp = sp_evaluate_element_list_param (space, msh, 'value', value, 'gradient', grad_param, 'hessian', hessian);

if (isempty (msh.elem_list))
  return
end

if (hessian)
  [Jinv, Jinv2] = geopdes_inv_der2__ (msh.geo_map_jac, msh.geo_map_der2);
  Jinv  = reshape (Jinv, [msh.ndim, msh.rdim, msh.nqn, 1, msh.nel]);
  JinvT = permute (Jinv, [2 1 3 4 5]);
  Jinv2 = reshape (Jinv2, [msh.ndim, msh.rdim, msh.rdim, msh.nqn, 1, msh.nel]);

  shg = reshape (sp.shape_function_gradients, [msh.ndim, 1, 1, msh.nqn, space.nsh_max, msh.nel]);
  shh = reshape (sp.shape_function_hessians, msh.ndim, msh.ndim, msh.nqn, space.nsh_max, msh.nel);
  shape_function_hessians = zeros (msh.ndim, msh.ndim, msh.nqn, space.nsh_max, msh.nel);

  shh_size = [1, 1, msh.nqn, space.nsh_max, msh.nel];
  for idim = 1:msh.rdim
    for jdim = 1:msh.rdim
      D2v_DF = sum (bsxfun(@times, shh, Jinv(:,idim,:,:,:)),1);
      DFT_D2v_DF = sum (bsxfun (@times, JinvT(jdim,:,:,:,:), D2v_DF), 2);
      Dv_D2F = sum (bsxfun (@times, shg, Jinv2(:,idim,jdim,:,:,:)), 1);

      shape_function_hessians(idim,jdim,:,:,:) = reshape (DFT_D2v_DF, shh_size) + ...
          reshape (Dv_D2F, shh_size);
    end
  end
  sp.shape_function_hessians = shape_function_hessians;
elseif (isfield (sp, 'shape_function_hessians'))
  sp = rmfield (sp, 'shape_function_hessians');
end

if (gradient)
  JinvT = geopdes_invT__ (msh.geo_map_jac);
  JinvT = reshape (JinvT, [msh.rdim, msh.ndim, msh.nqn, msh.nel]);
  sp.shape_function_gradients = geopdes_prod__ (JinvT, sp.shape_function_gradients);
elseif (isfield (sp, 'shape_function_gradients'))
  sp = rmfield (sp, 'shape_function_gradients');
end

end
