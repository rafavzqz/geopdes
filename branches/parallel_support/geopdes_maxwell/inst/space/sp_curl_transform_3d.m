% SP_CURL_TRANSFORM_3D: map a function space to the physical domain with a curl conserving transform.
%
%     sp = sp_curl_transform_3d (sp, msh)
%
% INPUTS:
%
%     sp:  space structure defined in the parametric domain (see sp_bspline_3d_param, sp_bspline_curl_transform_3d)
%     msh: structure containing the domain partition and the quadrature rule (see msh_push_forward_3d)
%
% OUTPUT:
%
%    sp: structure representing the discrete function space, with the following fields:
%
%        FIELD_NAME      (SIZE)                            DESCRIPTION
%        ndof            (scalar)                          total number of degrees of freedom
%        nsh_max         (scalar)                          maximum number of shape functions per element
%        nsh             (1 x msh.nel vector)              actual number of shape functions per each element
%        connectivity    (nsh_max x msh.nel vector)        indices of basis functions that do not vanish in each element
%        shape_functions (3 x msh.nqn x nsh_max x msh.nel) vector-valued basis functions evaluated at each quadrature node in each element
%        shape_function_curls
%                        (3 x msh.nqn x nsh_max x msh.nel) basis functions curls evaluated at each quadrature node in each element
%        boundary        (1 x 6 struct array)              struct array representing the space of tangential traces of basis functions on each edge
%        spfun           (function handle)                 function to evaluate an element of the discrete function space, given the Fourier coefficients and a set of points in the parametric space
%
%   For more details, see the documentation
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2010 Rafael Vazquez
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
  
function sp = sp_curl_transform_3d (sp, msh)

  sp = do_sp_curl_transform_3d__ (sp, msh);

  if (isfield (msh, 'boundary'))
    for iside = 1:numel(msh.boundary)
      sp.boundary(iside) = do_sp_curl_transform_3d__ (sp.boundary(iside), msh.boundary(iside));
    end
  end

end


function sp = do_sp_curl_transform_3d__ (sp, msh)

  [JinvT, jacdet] = geopdes_invT__ (msh.geo_map_jac);
  sp.shape_functions = geopdes_prod__ (JinvT, sp.shape_functions);

  if (isfield (sp, 'shape_function_curls'))
    DFcurl = geopdes_prod__ (msh.geo_map_jac, sp.shape_function_curls);
    for ii=1:sp.nsh_max
      DFaux = DFcurl(1,:,ii,:);
      sp.shape_function_curls(1,:,ii,:) = ...
          reshape(DFaux(:)./jacdet(:), 1, msh.nqn, 1, msh.nel);
      DFaux = DFcurl(2,:,ii,:);
      sp.shape_function_curls(2,:,ii,:) = ...
          reshape(DFaux(:)./jacdet(:), 1, msh.nqn, 1, msh.nel);
      DFaux = DFcurl(3,:,ii,:);
      sp.shape_function_curls(3,:,ii,:) = ...
          reshape(DFaux(:)./jacdet(:), 1, msh.nqn, 1, msh.nel);
    end
    clear DFcurl DFaux
  end

  if (isfield (sp, 'shape_function_gradients'))
    error ('sp_curl_transform_3d: the computation of the gradient with the curl transform is not allowed yet')
  end

  if (isfield (sp, 'shape_function_divergence'))
    error ('sp_curl_transform_3d: the computation of the divergence with the curl transform is not allowed yet')
  end

  clear JinvT jacdet

end
