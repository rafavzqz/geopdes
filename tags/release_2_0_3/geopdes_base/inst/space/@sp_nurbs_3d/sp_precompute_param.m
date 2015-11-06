% SP_PRECOMPUTE_PARAM: precompute all the fields, as in the space structure of the technical report, before mapping to the physical domain.
%  This function is used in vectorial spaces, before applying the map.
%
%     space = sp_precompute_param (space, msh); computes all the fields
%     space = sp_precompute_param (space, msh, 'option');  only computes the selected fields
%
% INPUT:
%     
%    space: object representing the discrete function space (see sp_nurbs_2d).
%    'option', value: additional optional parameters, available options are:
%        nsh, connectivity, value (shape_functions), gradient (shape_function_gradients).
%     The value must be true or false. All the values are false by default.
%
% OUTPUT:
%
%    space: object containing the information of the input object, plus the 
%            fields of the old structure, that are listed below. If no option
%            is given all the fields are computed. If an option is given,
%            only the selected fields will be computed.
%
%    FIELD_NAME      (SIZE)                             DESCRIPTION
%    nsh             (1 x msh.nel vector)               actual number of shape functions per each element
%    connectivity    (nsh_max x msh.nel vector)         indices of basis functions that do not vanish in each element
%    shape_functions (msh.nqn x nsh_max x msh.nel)      basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%             (3 x msh.nqn x nsh_max x msh.nel)         basis function gradients evaluated at each quadrature node in each element
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

  if (isempty (varargin))
    nsh = true;
    connectivity = true;
    value = true;
    gradient = true;
  else
    if (~rem (length (varargin), 2) == 0)
      error ('sp_precompute: options must be passed in the [option, value] format');
    end
    nsh = false;
    connectivity = false;
    value = false;
    gradient = false;
    for ii=1:2:length(varargin)-1
      if (strcmpi (varargin{ii}, 'connectivity'))
        connectivity = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'nsh'))
        nsh = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'value'))
        value = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'gradient'))
        gradient = varargin{ii+1};
      else
        error ('sp_precompute_param: unknown option %s', varargin {ii});
      end
    end    
  end

  nelu = msh.nel_dir(1); nelv = msh.nel_dir(2); nelw = msh.nel_dir(3);
  nel = msh.nel;
  nqnu = msh.nqn_dir(1); nqnv = msh.nqn_dir(2); nqnw = msh.nqn_dir(3);
  nqn = msh.nqn;

  spu = sp.spu; spv = sp.spv; spw = sp.spw;
  nsh_max = sp.nsh_max;

  if (nsh)
    [NSHV, NSHU, NSHW] = meshgrid (spv.nsh, spu.nsh, spw.nsh);
    sp.nsh  = reshape (NSHU .* NSHV .* NSHW, 1, []);
  end

  conn_u = reshape (spu.connectivity, spu.nsh_max, 1, 1, nelu, 1, 1);
  conn_u = repmat  (conn_u, [1, spv.nsh_max, spw.nsh_max, 1, nelv, nelw]);
  conn_u = reshape (conn_u, [], nel);

  conn_v = reshape (spv.connectivity, 1, spv.nsh_max, 1, 1, nelv, 1);
  conn_v = repmat  (conn_v, [spu.nsh_max, 1, spw.nsh_max, nelu, 1, nelw]);
  conn_v = reshape (conn_v, [], nel);

  conn_w = reshape (spw.connectivity, 1, 1, spw.nsh_max, 1, 1, nelw);
  conn_w = repmat  (conn_w, [spu.nsh_max, spv.nsh_max, 1, nelu, nelv, 1]);
  conn_w = reshape (conn_w, [], nel);

  connectivity = zeros (nsh_max, nel);
  indices = (conn_u ~= 0) & (conn_v ~= 0) & (conn_w~=0);
    connectivity(indices) = ...
       sub2ind ([spu.ndof, spv.ndof, spw.ndof], conn_u(indices), conn_v(indices), conn_w(indices));
  clear conn_u conn_v conn_w

  if (connectivity)
    sp.connectivity = reshape (connectivity, nsh_max, nel);
  end

  if (value || gradient)
    shp_u = reshape (spu.shape_functions, nqnu, 1, 1, spu.nsh_max, 1, 1, nelu, 1, 1);
    shp_u = repmat  (shp_u, [1, nqnv, nqnw, 1, spv.nsh_max, spw.nsh_max, 1, nelv, nelw]);
    shp_u = reshape (shp_u, nqn, nsh_max, nel);

    shp_v = reshape (spv.shape_functions, 1, nqnv, 1, 1, spv.nsh_max, 1, 1, nelv, 1);
    shp_v = repmat  (shp_v, [nqnu, 1, nqnw, spu.nsh_max, 1, spw.nsh_max, nelu, 1, nelw]);
    shp_v = reshape (shp_v, nqn, nsh_max, nel);

    shp_w = reshape (spw.shape_functions, 1, 1, nqnw, 1, 1, spw.nsh_max, 1, 1, nelw);
    shp_w = repmat  (shp_w, [nqnu, nqnv, 1, spu.nsh_max, spv.nsh_max, 1, nelu, nelv, 1]);
    shp_w = reshape (shp_w, nqn, nsh_max, nel);

    W = sp.weights(connectivity);
    W = repmat (reshape (W, 1, sp.nsh_max, msh.nel), [msh.nqn, 1, 1]);
    shape_functions = W .* shp_u .* shp_v .* shp_w;
    D = repmat (reshape (sum (shape_functions, 2), msh.nqn, 1, msh.nel), [1, sp.nsh_max, 1]);
    shape_functions = shape_functions ./ D;

    if (value)
      sp.shape_functions = shape_functions;
    end
  end

  if (gradient)
    shg_u = reshape (spu.shape_function_gradients, nqnu, 1, 1, spu.nsh_max, 1, 1, nelu, 1, 1);
    shg_u = repmat  (shg_u, [1, nqnv, nqnw, 1, spv.nsh_max, spw.nsh_max, 1, nelv, nelw]);
    shg_u = reshape (shg_u, nqn, nsh_max, nel);
  
    shg_v = reshape (spv.shape_function_gradients, 1, nqnv, 1, 1, spv.nsh_max, 1, 1, nelv, 1);
    shg_v = repmat  (shg_v, [nqnu, 1, nqnw, spu.nsh_max, 1, spw.nsh_max, nelu, 1, nelw]);
    shg_v = reshape (shg_v, nqn, nsh_max, nel);

    shg_w = reshape (spw.shape_function_gradients, 1, 1, nqnw, 1, 1, spw.nsh_max, 1, 1, nelw);
    shg_w = repmat  (shg_w, [nqnu, nqnv, 1, spu.nsh_max, spv.nsh_max, 1, nelu, nelv, 1]);
    shg_w = reshape (shg_w, nqn, nsh_max, nel);

    Bu = W .* shg_u .* shp_v .* shp_w;
    Bv = W .* shp_u .* shg_v .* shp_w;
    Bw = W .* shp_u .* shp_v .* shg_w;

    clear shg_u shg_v shg_w

    Du = repmat (reshape (sum (Bu, 2), msh.nqn, 1, msh.nel), [1, sp.nsh_max, 1]);
    Dv = repmat (reshape (sum (Bv, 2), msh.nqn, 1, msh.nel), [1, sp.nsh_max, 1]);
    Dw = repmat (reshape (sum (Bw, 2), msh.nqn, 1, msh.nel), [1, sp.nsh_max, 1]);

    sp.shape_function_gradients(1,:,:,:) = (Bu - shape_functions .* Du)./D;
    sp.shape_function_gradients(2,:,:,:) = (Bv - shape_functions .* Dv)./D;
    sp.shape_function_gradients(3,:,:,:) = (Bw - shape_functions .* Dw)./D;
    sp.shape_function_gradients = reshape (sp.shape_function_gradients, [3, msh.nqn, sp.nsh_max, msh.nel]);

    clear Bu Bv Bw Du Dv Dw shape_functions
  end

end
