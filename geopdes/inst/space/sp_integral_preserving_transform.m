% SP_INTEGRAL_PRESERVING_TRANSFORM: apply the integral-preserving transform to the functions in the parametric domain
%
%     sp = sp_integral_preserving_transform (space, msh)
%
% INPUTS:
%     
%    space:   structure with the information in the parametric domain (see sp_scalar/sp_evaluate_col)
%    msh:     msh structure containing the information of the parametrization
%              in the points where basis functions have to be computed (see msh_cartesian/msh_evaluate_col)
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
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

function sp = sp_integral_preserving_transform (sp, msh, value)

  if (nargin < 3 && isfield (sp, 'shape_functions'))
    value = true;
  end

  if (value)
    jacdet = reshape (geopdes_det__ (msh.geo_map_jac), msh.nqn, 1, msh.nel);
    sp.shape_functions = bsxfun (@rdivide, sp.shape_functions, jacdet);
    sp.shape_functions = sp.shape_functions;
    if (isfield (msh, 'side_number'))
      sp.shape_functions = sp.shape_functions * (-1)^msh.side_number;
    end
  end
  
  if (isfield (sp, 'shape_function_gradients'))
    sp = rmfield (sp, 'shape_function_gradients');
  end
  if (isfield (sp, 'shape_function_hessians'))
    sp = rmfield (sp, 'shape_function_hessians');
  end

end