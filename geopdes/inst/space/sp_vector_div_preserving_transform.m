% SP_VECTOR_DIV_PRESERVING_TRANSFORM: apply the div-preserving transform to the functions in the parametric domain
%
%     sp = sp_vector_div_preserving_transform (space, msh, value, gradient, curl, divergence)
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

function sp = sp_vector_div_preserving_transform (sp, msh, value, gradient, curl, divergence)

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

  jacdet = reshape (geopdes_det__ (msh.geo_map_jac), msh.nqn, msh.nel);
  
  if (divergence  && isfield (sp, 'shape_function_divs'))
    jacdet = reshape (jacdet, msh.nqn, 1, msh.nel);
    sp.shape_function_divs = bsxfun (@rdivide, sp.shape_function_divs, jacdet);
  end

  if (gradient || curl || (divergence && ~isfield (sp, 'shape_function_divs')))
    shape_fun_grads = piola_transform_grad__ (sp, msh, sp.shape_functions, sp.shape_function_gradients);
  
    if (gradient)
      sp.shape_function_gradients = shape_fun_grads;
    end
    
    if (divergence && ~isfield (sp, 'shape_function_divs'))
      sp.shape_function_divs = zeros (msh.nqn, sp.nsh_max, msh.nel);
      for icomp = 1:sp.ncomp
        sp.shape_function_divs = sp.shape_function_divs + sp.shape_function_gradients(icomp,icomp,:,:,:);
      end
    end

    if (curl)
      if (msh.ndim ~= msh.rdim)
        error ('The computation of the curl has not been implemented for manifolds')
      end
      if (sp.ncomp == 2)
        sp.shape_function_curls = reshape (shape_fun_grads(2,1,:,:,:) - ...
	 	       shape_fun_grads(1,2,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
      elseif (sp.ncomp == 3)
        shape_fun_curls = zeros (sp.ncomp, msh.nqn, sp.nsh_max, msh.nel);
        for icomp = 1:sp.ncomp
          ind1 = mod(icomp,3) + 1;
          ind2 = mod(ind1, 3) + 1;
          shape_fun_curls(icomp,:,:,:) = reshape (shape_fun_grads(ind2,ind1,:,:,:) - ...
                     shape_fun_grads(ind1,ind2,:,:,:), 1, msh.nqn, sp.nsh_max, msh.nel);
        end
        sp.shape_function_curls = shape_fun_curls;
      end
    end
  end

% The transformation is applied to the value only at the end, to avoid
%  conflicts with the computation of the gradient
  if (value)
    sp.shape_functions = geopdes_prod__ (msh.geo_map_jac, sp.shape_functions);
    jacdet = reshape (jacdet, 1, msh.nqn, 1, msh.nel);
    sp.shape_functions = bsxfun (@rdivide, sp.shape_functions, jacdet);
  end
  
  if (~gradient && isfield (sp, 'shape_function_gradients'))
    sp = rmfield (sp, 'shape_function_gradients');
  end

  if (~value && isfield (sp, 'shape_functions'))
    sp = rmfield (sp, 'shape_functions');
  end
  
  if (isfield (sp, 'shape_function_hessians'))
    sp = rmfield (sp, 'shape_function_hessians');
  end
  
end



function shape_fun_grads = piola_transform_grad__ (sp, msh, shp, shg)
% PIOLA_TRANSFORM_GRAD__: compute the gradients for the div-conforming Piola transform.
%
%     shape_fun_grads = piola_transform_grad__ (space, msh, shp, shg)
%
% INPUT:
%     
%    space: structure representing the discrete function space (see sp_vector/sp_evaluate_col).
%    msh: structure containing the quadrature information (see msh_cartesian/msh_evaluate_col)
%    shp: basis functions evaluated in the parametric domain
%    shg: gradient of the basis functions, in the parametric domain
%
% OUTPUT:
%
%    shape_fun_grads: gradient of the basis functions in the physical domain
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

  shape_fun_grads = zeros (sp.ncomp, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);

  DF  = reshape (msh.geo_map_jac, msh.rdim, msh.ndim, 1, msh.nqn, 1, msh.nel);
  DF2 = reshape (msh.geo_map_der2, msh.rdim, msh.ndim, msh.ndim, msh.nqn, 1, msh.nel);
  shp = reshape (shp, 1, sp.ncomp_param, 1, msh.nqn, sp.nsh_max, msh.nel);
  shg = reshape (shg, 1, sp.ncomp_param, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);

  Jinv = geopdes_inv__ (msh.geo_map_jac);
  Jinv  = reshape (Jinv, [1, msh.ndim, msh.rdim, msh.nqn, 1, msh.nel]);
  JinvT = permute (Jinv, [3 2 1 4 5 6]);
  jacdet = reshape (geopdes_det__ (msh.geo_map_jac), msh.nqn, 1, msh.nel);


% Derivatives of |DF|, divided by |DF| to simplify
  Ddet = zeros (msh.rdim, msh.nqn, 1, msh.nel);
  DF2_JinvT = sum (sum (bsxfun (@times, DF2, JinvT), 1), 2);
  JinvT = reshape (JinvT, [1, msh.rdim, msh.ndim, msh.nqn, 1, msh.nel]);
  for idim = 1:msh.rdim
    aux = JinvT(:,idim,:,:,:,:);
    Ddet(idim,:,:,:) = sum (bsxfun (@times, aux, DF2_JinvT), 3);
  end
  clear aux DF2_JinvT

  for icomp = 1:sp.ncomp
    aux2 = sum (bsxfun (@times, DF(icomp,:,:,:,:,:), shg), 2);
    aux3 = reshape (sum (bsxfun (@times, DF(icomp,:,:,:,:,:), shp), 2), [1, msh.nqn, sp.nsh_max, msh.nel]);

    for jdim = 1:msh.rdim
      aux = sum (bsxfun (@times, DF2(icomp,:,:,:,:,:), JinvT(:,jdim,:,:,:,:)), 3);
      term1 = reshape (sum (bsxfun (@times, aux, shp), 2), [msh.nqn, sp.nsh_max, msh.nel]);

      term2 = reshape (sum (bsxfun (@times, aux2, JinvT(:,jdim,:,:,:,:)), 3), [msh.nqn, sp.nsh_max, msh.nel]);
      
      term3 = reshape (bsxfun (@times, -Ddet(jdim,:,:,:), aux3), [msh.nqn, sp.nsh_max, msh.nel]);

      shape_fun_grads(icomp,jdim,:,:,:) = bsxfun (@rdivide, term1 + term2 + term3, jacdet);
    end
  end

end
