% MSH_PRECOMPUTE: precompute all the fields, as in the msh structure of the technical report.
%
%     msh = msh_precompute (msh);
%     msh = msh_precompute (msh, 'option', value);
%
% INPUTS:
%     
%     msh: mesh object containing the quadrature information (see msh_2d)
%    'option', value: additional optional parameters, available options are:
%        quad_nodes, quad_weights, geo_map, geo_map_jac, geo_map_der2, jacdet.
%     The value must be true or false. All the values are false by default.
%
% OUTPUT:
%
%     msh: mesh object containing the information of the input object, plus the 
%            fields of the old structure, that are listed below. If no option
%            is given all the fields are computed. If an option is given,
%            only the selected fields will be computed.
%
%     FIELD_NAME    (SIZE)                  DESCRIPTION
%
%     quad_nodes    (2 x nqn x nel vector)  coordinates of the quadrature nodes in parametric space
%     quad_weights  (nqn x nel vector)      weights associated to the quadrature nodes
%     geo_map       (2 x nqn x nel vector)  physical coordinates of the quadrature nodes
%     geo_map_jac   (2 x 2 x nqn x nel)     Jacobian matrix of the map evaluated at the quadrature nodes
%     geo_map_der2  (2 x 2 x 2 x nqn x nel) Second order derivatives of the map evaluated at the quadrature nodes ()
%     jacdet        (nqn x nel)             determinant of the Jacobian evaluated in the quadrature points
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

function msh = msh_precompute (msh, varargin)

  if (nargin == 1)
    quad_nodes = true;
    quad_weights = true;
    geo_map = true;
    geo_map_jac = true;
    geo_map_der2 = true;
    jacdet = true;
  else
    if (~rem (length (varargin), 2) == 0)
      error ('msh_precompute: options must be passed in the [option, value] format');
    end
    quad_nodes = false;
    quad_weights = false;
    geo_map = false;
    geo_map_jac = false;
    geo_map_der2 = false;
    jacdet = false;
    for ii=1:2:length(varargin)-1
      if (strcmpi (varargin{ii}, 'quad_nodes'))
        quad_nodes = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'quad_weights'))
        quad_weights = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'geo_map'))
        geo_map = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'geo_map_jac'))
        geo_map_jac = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'geo_map_der2'))
        geo_map_der2 = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'jacdet'))
        jacdet = varargin{ii+1};
      else
        error ('msh_precompute: unknown option %s', varargin {ii});
      end
    end    
  end

  nelu = msh.nel_dir(1); nelv = msh.nel_dir(2);
  nel  = msh.nel;
  nqnu = msh.nqn_dir(1); nqnv = msh.nqn_dir(2);
  nqn  = msh.nqn;
  qn = msh.qn;

  if (quad_nodes || geo_map)
    quad_nodes_u = reshape (qn{1}, nqnu, 1, nelu, 1);
    quad_nodes_u = repmat  (quad_nodes_u, [1, nqnv, 1, nelv]);
    quad_nodes_u = reshape (quad_nodes_u, [], nel);

    quad_nodes_v = reshape (qn{2}, 1, nqnv, 1, nelv);
    quad_nodes_v = repmat  (quad_nodes_v, [nqnu, 1, nelu, 1]);
    quad_nodes_v = reshape (quad_nodes_v, [], nel);

    if (quad_nodes)
      msh.quad_nodes(1, :, :) = quad_nodes_u;
      msh.quad_nodes(2, :, :) = quad_nodes_v;
    end
    if (geo_map) % The map is applied below
      msh.geo_map(1, :, :) = quad_nodes_u;
      msh.geo_map(2, :, :) = quad_nodes_v;
    end
    clear quad_nodes_u quad_nodes_v

  end

  if (~isempty (msh.qw) && quad_weights)
    qwu = msh.qw{1};  qwv = msh.qw{2};
    quad_weights_u = reshape (qwu, nqnu, 1, nelu, 1);
    quad_weights_u = repmat  (quad_weights_u, [1, nqnv, 1, nelv]);
    quad_weights_u = reshape (quad_weights_u, [], nel);

    quad_weights_v = reshape (qwv, 1, nqnv, 1, nelv);
    quad_weights_v = repmat  (quad_weights_v, [nqnu, 1, nelu, 1]);
    quad_weights_v = reshape (quad_weights_v, [], nel);

    msh.quad_weights = quad_weights_u .* quad_weights_v;
    if (abs (sum (msh.quad_weights(:)) - 1) > 1e-10)
      warning ('msh_precompute: inconsistent quadrature formula')
    end
    clear quad_weights_u quad_weights_v
  end

  if (geo_map)
    F = feval (msh.map, {qn{1}(:)', qn{2}(:)'});
    F = reshape (F, [2, nqnu, nelu, nqnv, nelv]);
    F = permute (F, [1 2 4 3 5]);
    msh.geo_map = reshape (F, [2, nqn, nel]);
  end

  if (geo_map_jac || jacdet)
    jac = feval (msh.map_der, {qn{1}(:)', qn{2}(:)'});
    jac = reshape (jac, [2, 2, nqnu, nelu, nqnv, nelv]);
    jac = permute (jac, [1 2 3 5 4 6]);
    jac = reshape (jac, 2, 2, nqn, nel);
    if (geo_map_jac)
      msh.geo_map_jac = jac;
    end
    if (jacdet)
      msh.jacdet = abs (geopdes_det__ (jac));
      msh.jacdet = reshape (msh.jacdet, [nqn, nel]);
    end
  end

  if (~isempty (msh.map_der2) && geo_map_der2)
    hess = feval (msh.map_der2, {qn{1}(:)', qn{2}(:)'});
    hess = reshape (hess, [2, 2, 2, nqnu, nelu, nqnv, nelv]);
    hess = permute (hess, [1 2 3 4 6 5 7]);
    msh.geo_map_der2 = reshape (hess, 2, 2, 2, nqn, nel);
  end

end
