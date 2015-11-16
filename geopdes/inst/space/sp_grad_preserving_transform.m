% SP_GRAD_PRESERVING_TRANSFORM: apply the grad-preserving transform to the functions in the parametric domain
%
%     sp = sp_grad_preserving_transform (space, msh, [value, gradient, hessian, laplacian])
%
% INPUTS:
%     
%    space:   structure with the information in the parametric domain (see sp_scalar/sp_evaluate_col)
%    msh:     msh structure containing the information of the parametrization
%              in the points where basis functions have to be computed (see msh_cartesian/msh_evaluate_col)
%    value, gradient, hessian, laplacian: additional optional parameters, either true or false
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
%    sp: struct representing the discrete function space, with the following fields:
%              (see the article for a detailed description)
%
%    FIELD_NAME      (SIZE)                                 DESCRIPTION
%    ncomp           (scalar)                               number of components of the functions of the space (actually, 1)
%    ndof            (scalar)                               total number of degrees of freedom
%    ndof_dir        (ncomp x ndim matrix)                  for each component, number of degrees of freedom along each direction
%    nsh_max         (scalar)                               maximum number of shape functions per element
%    nsh             (1 x msh_col.nel vector)               actual number of shape functions per each element
%    connectivity    (nsh_max x msh_col.nel vector)         indices of basis functions that do not vanish in each element
%    shape_functions (msh_col.nqn x nsh_max x msh_col.nel)  basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%       (rdim x msh_col.nqn x nsh_max x msh_col.nel)        basis function gradients evaluated at each quadrature node in each element
%    shape_function_hessians
%       (rdim x rdim x msh_col.nqn x nsh_max x msh_col.nel) basis function hessians evaluated at each quadrature node in each element
%    shape_function_laplacians 
%       (msh_col.nqn x nsh_max x msh_col.nel)               basis functions laplacians evaluated at each quadrature node in each element
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

function sp = sp_grad_preserving_transform (sp, msh, value, gradient, hessian, laplacian)

  if (nargin < 3 || isempty (value))
    value = true;
  end
  if (nargin < 4 || isempty (gradient))
    gradient = false;
  end
  if (nargin < 5 || isempty (hessian))
    hessian = false;
  end
  if (nargin < 6 || isempty (laplacian))
    laplacian = false;
  end


  if (hessian || laplacian)
    [Jinv, Jinv2] = geopdes_inv_der2__ (msh.geo_map_jac, msh.geo_map_der2);
    Jinv  = reshape (Jinv, [msh.ndim, msh.rdim, msh.nqn, 1, msh.nel]);
    JinvT = permute (Jinv, [2 1 3 4 5]);
    Jinv2 = reshape (Jinv2, [msh.ndim, msh.rdim, msh.rdim, msh.nqn, 1, msh.nel]);

    shg = reshape (sp.shape_function_gradients, [msh.ndim, 1, 1, msh.nqn, sp.nsh_max, msh.nel]);
    shh = reshape (sp.shape_function_hessians, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);
    shape_function_hessians = zeros (msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);

    shh_size = [1, 1, msh.nqn, sp.nsh_max, msh.nel];
    for idim = 1:msh.rdim
      for jdim = 1:msh.rdim
        D2v_DF = sum (bsxfun(@times, shh, Jinv(:,idim,:,:,:)),1);
        DFT_D2v_DF = sum (bsxfun (@times, JinvT(jdim,:,:,:,:), D2v_DF), 2);
        Dv_D2F = sum (bsxfun (@times, shg, Jinv2(:,idim,jdim,:,:,:)), 1);

        shape_function_hessians(idim,jdim,:,:,:) = reshape (DFT_D2v_DF, shh_size) + ...
            reshape (Dv_D2F, shh_size);
      end
    end
      
    if (hessian)
      sp.shape_function_hessians = shape_function_hessians;
    end
      
    if (laplacian)
      sp.shape_function_laplacians = zeros (msh.nqn, sp.nsh_max, msh.nel);
      for idim = 1:msh.rdim
        sp.shape_function_laplacians = sp.shape_function_laplacians + ...
            reshape (shape_function_hessians(idim,idim,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
      end
    end
  end
  if (~hessian && isfield (sp, 'shape_function_hessians'))
      sp = rmfield (sp, 'shape_function_hessians');
  end
  if (~laplacian && isfield (sp, 'shape_function_laplacians'))
      sp = rmfield (sp, 'shape_function_laplacians');
  end

  if (gradient)
    JinvT = geopdes_invT__ (msh.geo_map_jac);
    JinvT = reshape (JinvT, [msh.rdim, msh.ndim, msh.nqn, msh.nel]);
    sp.shape_function_gradients = geopdes_prod__ (JinvT, sp.shape_function_gradients);
  elseif (isfield (sp, 'shape_function_gradients'))
    sp = rmfield (sp, 'shape_function_gradients');
  end

  if (~value && isfield (sp, 'shape_functions'))
    sp = rmfield (sp, 'shape_functions');
  end
  
end