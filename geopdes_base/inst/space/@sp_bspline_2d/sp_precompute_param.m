% SP_PRECOMPUTE_PARAM: precompute all the fields, as in the space structure of the technical report, before mapping to the physical domain.
%  This function is used in vectorial spaces, before applying the map.
%
%     space = sp_precompute_param (space, msh)
%     space = sp_precompute_param (space, msh, 'option')
%
% INPUT:
%     
%    space: object representing the discrete function space (see sp_bspline_2d).
%
% OUTPUT:
%
%    space: object representing the discrete function space in the parametric domain, plus the following fields (or some of them):
%
%    FIELD_NAME      (SIZE)                                 DESCRIPTION
%    connectivity    (nsh_max x msh_col.nel vector)         indices of basis functions that do not vanish in each element
%    shape_functions (msh_col.nqn x nsh_max x msh_col.nel)  basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%             (2 x msh_col.nqn x nsh_max x msh_col.nel)     basis function gradients evaluated at each quadrature node in each element
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011 Rafael Vazquez
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

function sp = sp_precompute_param (sp, msh, varargin)

  if (nargin == 2)
    connectivity = true;
    value = true;
    gradient = true;
  else
    connectivity = false;
    value = false;
    gradient = false;
    for ii=1:length(varargin)
      if (strcmpi (varargin {ii}, 'connectivity'))
        connectivity = true;
      elseif (strcmpi (varargin {ii}, 'value'))
        value = true;
      elseif (strcmpi (varargin {ii}, 'gradient'))
        gradient = true;
      else
        error ('sp_precompute_param: unknown option %s', varargin {ii});
      end
    end    
  end

  nelu = msh.nel_dir(1); nelv = msh.nel_dir(2);
  nel = msh.nel;
  nqnu = msh.nqn_dir(1); nqnv = msh.nqn_dir(2);
  nqn = msh.nqn;

  spu = sp.spu; spv = sp.spv;
  nsh_max = sp.nsh_max;

  if (connectivity)
    conn_u = reshape (spu.connectivity, spu.nsh_max, 1, nelu, 1);
    conn_u = repmat  (conn_u, [1, spv.nsh_max, 1, nelv]);
    conn_u = reshape (conn_u, [], nel);

    conn_v = reshape (spv.connectivity, 1, spv.nsh_max, 1, nelv);
    conn_v = repmat  (conn_v, [spu.nsh_max, 1, nelu, 1]);
    conn_v = reshape (conn_v, [], nel);

    connectivity = zeros (nsh_max, nel);
    indices = (conn_u ~= 0) & (conn_v ~= 0);
    connectivity(indices) = ...
        sub2ind ([spu.ndof, spv.ndof], conn_u(indices), conn_v(indices));
    sp.connectivity = reshape (connectivity, nsh_max, nel);

    clear conn_u conn_v connectivity
  end

  if (value || gradient)
    shp_u = reshape (spu.shape_functions, nqnu, 1, spu.nsh_max, 1, nelu, 1);
    shp_u = repmat  (shp_u, [1, nqnv, 1, spv.nsh_max, 1, nelv]);
    shp_u = reshape (shp_u, nqn, nsh_max, nel);

    shp_v = reshape (spv.shape_functions, 1, nqnv, 1, spv.nsh_max, 1, nelv);
    shp_v = repmat  (shp_v, [nqnu, 1, spu.nsh_max, 1, nelu, 1]);
    shp_v = reshape (shp_v, nqn, nsh_max, nel);

    if (value)
      sp.shape_functions = shp_u .* shp_v;
    end
  end

  if (gradient)
    shg_u = reshape (spu.shape_function_gradients, nqnu, 1, spu.nsh_max, 1, nelu, 1);
    shg_u = repmat  (shg_u, [1, nqnv, 1, spv.nsh_max, 1, nelv]);
    shg_u = reshape (shg_u, nqn, nsh_max, nel);
  
    shg_v = reshape (spv.shape_function_gradients, 1, nqnv, 1, spv.nsh_max, 1, nelv);
    shg_v = repmat  (shg_v, [nqnu, 1, spu.nsh_max, 1, nelu, 1]);
    shg_v = reshape (shg_v, nqn, nsh_max, nel);

    sp.shape_function_gradients(1,:,:,:) = shg_u .* shp_v ;
    sp.shape_function_gradients(2,:,:,:) = shp_u .* shg_v ;
    sp.shape_function_gradients = reshape (sp.shape_function_gradients, [2, msh.nqn, sp.nsh_max, msh.nel]);

    clear shg_u shg_v
  end

end
