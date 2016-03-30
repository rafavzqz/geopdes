% SP_VECTOR_CURL_PRESERVING_TRANSFORM: apply the curl-preserving transform to the functions in the parametric domain
%
%     sp = sp_vector_curl_preserving_transform (space, msh, value, curl)
%
% INPUTS:
%     
%    space:   structure with the information in the parametric domain (see sp_vector/sp_evaluate_col_param)
%    msh:     msh structure containing the information of the parametrization
%              in the points where basis functions have to be computed (see msh_cartesian/msh_evaluate_col)
%    value, curl: additional optional parameters, either true or false
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
%            curl       |      false      |  compute shape_function_curls
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

function sp = sp_vector_curl_preserving_transform (sp, msh, value, curl)

  if (nargin < 3)
    value = true;
  end
  if (nargin < 4)
    curl = false;
  end

  [JinvT, jacdet] = geopdes_invT__ (msh.geo_map_jac);

  if (value)
    JinvT = reshape (JinvT, [msh.rdim, msh.ndim, msh.nqn, msh.nel]);
    sp.shape_functions = geopdes_prod__ (JinvT, sp.shape_functions);
  elseif (isfield (sp, 'shape_functions'))
    sp = rmfield (sp, 'shape_functions');
  end

  if (curl)
    if (sp.ncomp_param == 2)
      jacdet = reshape (jacdet, msh.nqn, 1, msh.nel);
      sp.shape_function_curls = bsxfun (@rdivide, sp.shape_function_curls, jacdet);
  
    elseif (sp.ncomp_param == 3)
      jacdet = reshape (jacdet, 1, msh.nqn, 1, msh.nel);
      sp.shape_function_curls = geopdes_prod__ (msh.geo_map_jac, sp.shape_function_curls);
      sp.shape_function_curls = bsxfun (@rdivide, sp.shape_function_curls, jacdet);
    end
  elseif (isfield (sp, 'shape_function_curls'))
    sp = rmfield (sp, 'shape_function_curls');
  end

  
  if (isfield (sp, 'shape_function_divs'))
    sp = rmfield (sp, 'shape_function_divs');
  end
    
  if (isfield (sp, 'shape_function_gradients'))
    sp = rmfield (sp, 'shape_function_gradients');
  end

  if (isfield (sp, 'shape_function_hessians'))
    sp = rmfield (sp, 'shape_function_hessians');
  end
  
end