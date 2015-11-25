% SP_VECTOR_GRAD_PRESERVING_TRANSFORM: apply the grad-preserving transform to the functions in the parametric domain
%
%     sp = sp_vector_grad_preserving_transform (space, msh, value, gradient, curl, divergence)
%
% INPUTS:
%     
%    space:   structure with the information in the parametric domain (see sp_vector/sp_evaluate_col_param)
%    msh:     msh structure containing the information of the parametrization
%              in the points where basis functions have to be computed (see msh_cartesian/msh_evaluate_col)
%    value, gradient, curl, divergence: additional optional parameters, either true or false
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
%            gradient   |      false      |  compute shape_function_gradients
%            curl       |      false      |  compute shape_function_curls
%            divergence |      false      |  compute shape_function_divs
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
%    nsh             (1 x msh_col.nel vector)                   actual number of shape functions per each element
%    connectivity    (nsh_max x msh_col.nel vector)             indices of basis functions that do not vanish in each element
%    shape_functions (ncomp x msh_col.nqn x nsh_max x msh_col.nel)  basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%       (ncomp x rdim x msh_col.nqn x nsh_max x msh_col.nel)    basis function gradients evaluated at each quadrature node in each element
%    shape_function_divs (msh_col.nqn x nsh_max x msh_col.nel)  basis function divergence evaluated at each quadrature node in each element
%    shape_function_curls 
%         2D:  (msh_col.nqn x nsh_max x msh_col.nel)            basis function curl evaluated at each quadrature node in each element
%         3D:  (3 x msh_col.nqn x nsh_max x msh_col.nel)        
%
% Copyright (C) 2015 Rafael Vazquez
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

function sp = sp_vector_grad_preserving_transform (sp, msh, value, gradient, curl, divergence)

  if (nargin < 3)
    value = true;
  end
  if (nargin < 4)
    gradient = false;
  end
  if (nargin < 5)
    curl = false;
  end
  if (nargin < 6)
    divergence = false;
  end

  if (gradient || curl || divergence)
    JinvT = geopdes_invT__ (msh.geo_map_jac);
    JinvT = reshape (JinvT, [msh.rdim, msh.ndim, msh.nqn, msh.nel]);
    grads = zeros (sp.ncomp, msh.rdim, msh.nqn, sp.nsh_max, msh.nel);
    for icomp = 1:sp.ncomp
      grads(icomp,:,:,:,:) = geopdes_prod__ (JinvT, ...
          reshape (sp.shape_function_gradients(icomp,:,:,:,:), msh.ndim, msh.nqn, sp.nsh_max, msh.nel));
    end
    sp.shape_function_gradients = grads;
    clear grads
    
    if (divergence)
      if (sp.ncomp == msh.rdim)
        sp.shape_function_divs = zeros (msh.nqn, sp.nsh_max, msh.nel);
        for icomp = 1:sp.ncomp
          sp.shape_function_divs = sp.shape_function_divs + ...
            reshape (sp.shape_function_gradients(icomp,icomp,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
        end
      else
        error('The divergence is not implemented in the case that rdim != ncomp')
      end
    elseif (isfield (sp, 'shape_function_divs'))
      sp = rmfield (sp, 'shape_function_divs');
    end

    if (curl)
      if (sp.ncomp == 2 && msh.rdim == 2)
        sp.shape_function_curls = reshape (sp.shape_function_gradients(2,1,:,:,:) - ...
			 	       sp.shape_function_gradients(1,2,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
      elseif (sp.ncomp == 3 && msh.rdim == 3)
        sp.shape_function_curls = zeros (sp.ncomp, msh.nqn, sp.nsh_max, msh.nel);
        for icomp = 1:sp.ncomp
          ind1 = mod(icomp,3) + 1;
          ind2 = mod(ind1, 3) + 1;
          sp.shape_function_curls(icomp,:,:,:) = reshape (sp.shape_function_gradients(ind2,ind1,:,:,:) - ...
                     sp.shape_function_gradients(ind1,ind2,:,:,:), 1, msh.nqn, sp.nsh_max, msh.nel);
        end
      else
        error('The curl is not implemented in the case that rdim != ncomp')
      end
    elseif (isfield (sp, 'shape_function_curls'))
      sp = rmfield (sp, 'shape_function_curls');
    end
    
    if (~gradient)
      sp = rmfield (sp, 'shape_function_gradients');
    end
  end

  if (~value && isfield (sp, 'shape_functions'))
    sp = rmfield (sp, 'shape_functions');
  end
  
end