% SP_VECTOR_GRAD_PRESERVING_TRANSFORM: apply the grad-preserving transform to the functions in the parametric domain
%
%     sp = sp_vector_grad_preserving_transform (space, msh, [value, gradient, curl, divergence, hessian])
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
%            hessian    |      false      |  compute shape_function_hessians
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

function sp = sp_vector_grad_preserving_transform (sp, msh, value, gradient, curl, divergence, hessian)

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
  if (nargin < 7)
    hessian = false;
  end

  if (hessian)
    [Jinv, Jinv2] = geopdes_inv_der2__ (msh.geo_map_jac, msh.geo_map_der2);
    Jinv  = reshape (Jinv, [msh.ndim, msh.rdim, msh.nqn, 1, msh.nel]);
    JinvT = permute (Jinv, [2 1 3 4 5]);
    Jinv2 = reshape (Jinv2, [msh.ndim, msh.rdim, msh.rdim, msh.nqn, 1, msh.nel]);

    shape_function_hessians = zeros (sp.ncomp, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);

    shh_size = [1, 1, 1, msh.nqn, sp.nsh_max, msh.nel];

    for icomp = 1:sp.ncomp
      shg = reshape (sp.shape_function_gradients(icomp,:,:,:,:), ...
        [msh.ndim, 1, 1, msh.nqn, sp.nsh_max, msh.nel]);
      shh = reshape (sp.shape_function_hessians(icomp,:,:,:,:,:), ...
         msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);
      for idim = 1:msh.rdim
        for jdim = 1:msh.rdim
          D2v_DF = sum (bsxfun(@times, shh, Jinv(:,idim,:,:,:)),1);
          DFT_D2v_DF = sum (bsxfun (@times, JinvT(jdim,:,:,:,:), D2v_DF), 2);
          Dv_D2F = sum (bsxfun (@times, shg, Jinv2(:,idim,jdim,:,:,:)), 1);

          shape_function_hessians(icomp, idim,jdim,:,:,:) = ...
            reshape (DFT_D2v_DF, shh_size) + reshape (Dv_D2F, shh_size);
        end
      end
    end

    if (hessian)
      sp.shape_function_hessians = shape_function_hessians;
    end
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
  end
  if (~gradient && isfield (sp, 'shape_function_gradients'))
    sp = rmfield (sp, 'shape_function_gradients');
  end

  if (~value && isfield (sp, 'shape_functions'))
    sp = rmfield (sp, 'shape_functions');
  end
  
end
