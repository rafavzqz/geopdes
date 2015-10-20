% MSH_PRECOMPUTE: precompute all the fields, as in the msh structure of the technical report.
%
%     msh = msh_precompute (msh);
%     msh = msh_precompute (msh, 'option', value);
%
% INPUTS:
%     
%     msh: mesh object containing the quadrature information (see msh_cartesian)
%    'option', value: additional optional parameters, available options are:
%        quad_nodes, quad_weights, geo_map, geo_map_jac, jacdet.
%     The value must be true or false. All the values are false by default.
%
% OUTPUT:
%
%     msh: mesh object containing the information of the input object, plus the 
%            fields of the old structure, that are listed below. If no option
%            is given all the fields are computed. If an option is given,
%            only the selected fields will be computed.
%
%     FIELD_NAME    (SIZE)                     DESCRIPTION
%
%     ndim          (scalar)                   dimension of the parametric space
%     rdim          (scalar)                   dimension of the physical space
%     nel           (scalar)                   number of elements in the mesh
%     nel_dir       (1 x ndim vector)          number of elements in each parametric direction for the entire mesh
%     nqn           (scalar)                   number of quadrature nodes per element
%     nqn_dir       (1 x ndim vector)          number of quadrature nodes per element in each parametric direction for the entire mesh
%     quad_nodes    (ndim x nqn x nel vector)  coordinates of the quadrature nodes in parametric space
%     quad_weights  (nqn x nel vector)         weights associated to the quadrature nodes
%     geo_map       (rdim x nqn x nel vector)  physical coordinates of the quadrature nodes
%     geo_map_jac   (rdim x ndim x nqn x nel)  Jacobian matrix of the map evaluated at the quadrature nodes
%     jacdet        (nqn x nel)                element of length, area, volume (if rdim = ndim, determinant of the Jacobian)
%     geo_map_der2  (rdim x ndim x ndim x nqn x nel]) Hessian matrix of the map evaluated at the quadrature nodes
%     normal        (rdim x ndim x nqn x nel]) for 3D surfaces, the exterior normal to the surface
%     element_size  (1 x nel)                  the size of the element in the physical domain
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011 Rafael Vazquez
% Copyright (C) 2014 Elena Bulgarello, Carlo de Falco, Sara Frizziero
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

function msh = msh_precompute (msh, varargin)

  if (nargin == 1)
    quad_nodes = true;
    quad_weights = true;
    geo_map = true;
    geo_map_jac = true;
    geo_map_der2 = true;
    jacdet = true;
    normal = true;
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
    normal = false;
    for ii=1:2:length(varargin)-1
      if (strcmpi (varargin{ii}, 'quad_nodes'))
        quad_nodes = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'quad_weights'))
        quad_weights = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'geo_map'))
        geo_map = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'geo_map_jac'))
        geo_map_jac = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'jacdet'))
        jacdet = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'geo_map_der2'))
        geo_map_der2 = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'normal'))
        normal = varargin{ii+1};
      else
        error ('msh_precompute: unknown option %s', varargin {ii});
      end
    end
    
  end

  msh = struct (msh);

  if (quad_nodes)
    for idim = 1:msh.ndim
      qsize = ones (1, msh.ndim*2);
%       qsize(2*idim-1:2*idim) = [msh.nqn_dir(idim), msh.nel_dir(idim)];
      qsize([idim, msh.ndim+idim]) = [msh.nqn_dir(idim), msh.nel_dir(idim)];
      qrep = [msh.nqn_dir(1:msh.ndim), msh.nel_dir(1:msh.ndim)];
      qrep([idim, msh.ndim+idim]) = 1;
      quad_nodes = reshape (msh.qn{idim}, qsize);
      quad_nodes = repmat (quad_nodes, qrep);
      msh.quad_nodes(idim,:,:) = reshape (quad_nodes, msh.nqn, msh.nel);
    end
    clear quad_nodes
  end

  if (~isempty (msh.qw) && quad_weights)
    qw = 1;
    for idim = 1:msh.ndim
      qw = kron (msh.qw{idim}, qw);
    end
    msh.quad_weights = qw;
    clear qw
  end

% Auxiliary vector sizes, to use with reshape and permute
  reorder = @(x) x(:)';
  psize = reorder ([msh.nqn_dir; msh.nel_dir]);
  vorder = [1:2:msh.ndim*2, 2:2:msh.ndim*2]; % [1 3 5 2 4 6], for ndim = 3
  if (geo_map)
    F = feval (msh.map, cellfun (reorder, msh.qn, 'UniformOutput', false));
    F = reshape (F, [msh.rdim, psize]);
    F = permute (F, [1, vorder+1]);
    msh.geo_map = reshape (F, [msh.rdim, msh.nqn, msh.nel]);
  end

  if (geo_map_jac || jacdet)
    jac = feval (msh.map_der, cellfun (reorder, msh.qn, 'UniformOutput', false));
    jac = reshape (jac, [msh.rdim, msh.ndim, psize]);
    jac = permute (jac, [1 2 vorder+2]);
    jac = reshape (jac, msh.rdim, msh.ndim, msh.nqn, msh.nel);

    if (geo_map_jac)
      msh.geo_map_jac = jac;
    end
    if (jacdet)
      msh.jacdet = abs (geopdes_det__ (jac));
      msh.jacdet = reshape (msh.jacdet, [msh.nqn, msh.nel]);
    end
  end

  if (~isempty (msh.map_der2) && geo_map_der2)
    hess = feval (msh.map_der2, cellfun (reorder, msh.qn, 'UniformOutput', false));
    hess = reshape (hess, [msh.rdim, msh.ndim, msh.ndim, psize]);
    hess = permute (hess, [1 2 3 vorder+3]);
    msh.geo_map_der2 = reshape (hess, [msh.rdim, msh.ndim, msh.ndim, msh.nqn, msh.nel]);
  end

  if (msh.ndim == 2 && msh.rdim == 3 && normal)
    normals = reshape (geopdes_cross__ (msh.geo_map_jac(:,1,:,:), ...
                              msh.geo_map_jac(:,2,:,:)), msh.rdim, msh.nqn, msh.nel);
    norms = reshape (geopdes_norm__ (normals), [1, msh.nqn, msh.nel]);
    msh.normal = bsxfun (@rdivide, normals, norms);
  end

%   if (any (ismember(fieldnames(msh), 'quad_weights')) && any (ismember(fieldnames(msh), 'jacdet')))
%     msh.element_size = (sum (msh.quad_weights .* ...
%                              abs (msh.jacdet), 1)).^(1/msh.ndim);
%   end

end
