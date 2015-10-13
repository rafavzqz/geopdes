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

function sp = sp_precompute (sp, msh, varargin)

  if (isempty (varargin))
    gradient = true;
    hessian = true;
    if (msh.ndim ~= msh.rdim)
      hessian = false;
    end
  else
    if (~rem (length (varargin), 2) == 0)
      error ('sp_precompute: options must be passed in the [option, value] format');
    end

    gradient = false;
    hessian = false;
    for ii=1:2:length(varargin)-1
      if (strcmpi (varargin{ii}, 'gradient'))
        gradient = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'hessian'))
        hessian = varargin{ii+1};
      end
    end    
  end

  if (isempty (varargin))
    sp = sp_precompute_param (sp, msh);
  else
    sp = sp_precompute_param (sp, msh, varargin{:});
  end

  switch (lower (space.transform))
    case {'grad-preserving'}
      sp = sp_grad_preserving_transform (sp, msh, value, gradient, hessian, laplacian);
    case {'integral-preserving'}
      sp = sp_integral_preserving_transform (sp, msh, value);
      if (gradient || hessian || laplacian)
        warning ('Derivatives not implemented for integral-preserving transformation')
      end
  end
  
  if (hessian)
    if (msh.rdim == msh.ndim)
      if (isempty (msh.geo_map_jac))
        msh = msh_precompute (msh, 'geo_map_jac', true);
      end
      if (isempty (msh.geo_map_der2))
        msh = msh_precompute (msh, 'geo_map_der2', true);
      end

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
      sp.shape_function_hessians = shape_function_hessians;
%     else
%       warning ('For manifolds, the second derivatives are not implemented yet')
    end
  elseif (isfield (sp, 'shape_function_hessians'))
    sp = rmfield (sp, 'shape_function_hessians');
  end
  
  if (gradient)
    if (isempty (msh.geo_map_jac))
      msh = msh_precompute (msh, 'geo_map_jac', true);
    end
    JinvT = geopdes_invT__ (msh.geo_map_jac);
    JinvT = reshape (JinvT, [msh.rdim, msh.ndim, msh.nqn, msh.nel]);
    sp.shape_function_gradients = geopdes_prod__ (JinvT, sp.shape_function_gradients);
  elseif (isfield (sp, 'shape_function_gradients'))
    sp = rmfield (sp, 'shape_function_gradients');
  end

end
