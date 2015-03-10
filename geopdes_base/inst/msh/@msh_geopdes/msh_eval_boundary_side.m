% MSH_EVAL_BOUNDARY_SIDE: evaluate the parameterization in one boundary side of the domain.
%
%     msh_side = msh_eval_boundary_side (msh, iside);
%
% INPUTS:
%     
%    msh:   mesh object (see msh_3d)
%    iside: number of the boundary side to compute, from 1 to 6 (see the file geo_specs for documentation about face numbering)
%
% OUTPUT:
%
%     msh_side: structure that contains the following fields
%
%     FIELD_NAME    (SIZE)                  DESCRIPTION
%     side_number   (scalar)                number of the side
%     nel           (scalar)                number of elements of the boundary side
%     nel_dir       (1 x 2 vector)          number of elements in each parametric direction
%     nqn           (scalar)                number of quadrature nodes per element
%     nqn_dir       (1 x 2 vector)          number of quadrature nodes per element in each parametric direction
%     quad_nodes    (3 x nqn x nel vector)  coordinates of the quadrature nodes in parametric domain
%     quad_weights  (nqn x nel vector)      weights associated to the quadrature nodes
%     geo_map       (3 x nqn x nel vector)  physical coordinates of the quadrature nodes
%     geo_map_jac   (3 x 3 x nqn x nel)     Jacobian matrix of the map evaluated at the quadrature nodes
%     jacdet        (nqn x nel)             element of length, area, volume (if rdim = ndim, determinant of the Jacobian)
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2014, 2015 Rafael Vazquez
% Copyright (C) 2014 Adriano Cortes
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

function msh_side = msh_eval_boundary_side (msh, iside)

  ind = setdiff (1:msh.ndim, ceil(iside/2)); % ind = [2 3; 2 3; 1 3; 1 3; 1 2; 1 2] 
  ind2 = floor ((iside+1)/2); % ind2 = [1 1 2 2 3 3];

  msh_side = msh.boundary(iside);

  qw = 1;
  for idim = 1:msh_side.ndim
    qw = kron (msh_side.qw{idim}, qw);
  end
  msh_side.quad_weights = qw;

  msh_side.quad_nodes = zeros ([msh.ndim, msh_side.nqn, msh_side.nel]);
  for idim = 1:msh_side.ndim
    qsize = ones (1, msh_side.ndim*2);
    qsize(2*idim-1:2*idim) = [msh_side.nqn_dir(idim), msh_side.nel_dir(idim)];
    qrep = [msh_side.nqn_dir(1:msh_side.ndim), msh_side.nel_dir(1:msh_side.ndim)];
    qrep([idim, msh_side.ndim+idim]) = 1;
    quad_nodes = reshape (msh_side.qn{idim}, qsize);
    quad_nodes = repmat (quad_nodes, qrep);
    msh_side.quad_nodes(ind(idim),:,:) = reshape (quad_nodes, msh_side.nqn, msh_side.nel);
  end
  clear quad_nodes

  qn(ind) = msh_side.qn;
  if (mod (iside, 2) == 0)
    msh_side.quad_nodes(ind2,:,:) = msh.breaks{ind2}(end);
    qn{ind2} = msh.breaks{ind2}(end);
  else
    msh_side.quad_nodes(ind2,:,:) = msh.breaks{ind2}(1);
    qn{ind2} = msh.breaks{ind2}(1);
  end

  
% Auxiliary vector sizes, to use with reshape and permute
  reorder = @(x) x(:)';
  psize = reorder ([msh_side.nqn_dir; msh_side.nel_dir]);
  vorder = [1:2:msh_side.ndim*2, 2:2:msh_side.ndim*2]; % [1 3 5 2 4 6], for ndim = 3
  
  F = feval (msh.map, cellfun (reorder, qn, 'UniformOutput', false));
  F = reshape (F, [msh.rdim, psize]);
  F = permute (F, [1, vorder+1]);
  jac = feval (msh.map_der, cellfun (reorder, qn, 'UniformOutput', false));
  jac = reshape (jac, [msh.rdim, msh.ndim, psize]);
  jac = permute (jac, [1 2 vorder+2]);

  msh_side.geo_map = reshape (F, msh_side.rdim, msh_side.nqn, msh_side.nel);
  msh_side.geo_map_jac = reshape(jac, msh.rdim, msh.ndim, msh_side.nqn, msh_side.nel);

%   jacdet = geopdes_norm__ (cross ...
%                           (squeeze (msh_side.geo_map_jac(:,ind(1),:,:)), ...
%                            squeeze (msh_side.geo_map_jac(:,ind(2),:,:))));
%   msh_side.jacdet = reshape (jacdet, msh_side.nqn, msh_side.nel);
  msh_side.jacdet = geopdes_det__ (msh_side.geo_map_jac(:,ind,:,:));

  if (msh.ndim == 2 || msh.ndim == 3)
% Compute the normal vector
    JinvT = geopdes_invT__ (msh_side.geo_map_jac);
    JinvT = reshape (JinvT, [msh.rdim, msh.ndim, msh_side.nqn, msh_side.nel]);

    normal = zeros (msh.ndim, msh_side.nqn, 1, msh_side.nel);
    normal(ind2,:) = (-1)^iside;
    normal = geopdes_prod__ (JinvT, normal);
    normal = reshape (normal, [msh.rdim, msh_side.nqn, msh_side.nel]);
  %Now normalize
    norms = reshape (geopdes_norm__ (normal), [1, msh_side.nqn, msh_side.nel]);
    msh_side.normal = bsxfun (@rdivide, normal, norms);

    if (msh.ndim == msh.rdim)
% Compute the characteristic length. This part should be improved
      charlen_param(ind) = cellfun (@(x) diff(x) * 0.5, msh.breaks(ind), 'UniformOutput', false);
      if (mod (iside,2) == 0)
        charlen_param{ind2} = 0.5*(msh.breaks{ind2}(end) - msh.breaks{ind2}(end-1));
      else
        charlen_param{ind2} = 0.5*(msh.breaks{ind2}(2) - msh.breaks{ind2}(1));
      end
      [charlen_param{1:msh.ndim}] = ndgrid (charlen_param{:});

      charlen = zeros (1, msh.ndim, 1, msh_side.nel);
      for idim = 1:msh.ndim
        charlen(:,idim,:,:) = reshape (charlen_param{idim}, 1, msh_side.nel);
      end

      geo_map_jac = bsxfun (@times, msh_side.geo_map_jac, charlen);
%       for iel = 1:msh_side.nel
%         geo_map_jac(:,1,:,iel) = xi_span_charlen(iel).*geo_map_jac(:,1,:,iel);
%         geo_map_jac(:,2,:,iel) = eta_span_charlen(iel).*geo_map_jac(:,2,:,iel);
%       end

      [Jinv, ~] = geopdes_inv__ (geo_map_jac);
      Jinv = reshape (Jinv, [msh.rdim, msh.ndim, msh_side.nqn, msh_side.nel]);
      normal = reshape (normal, [msh.rdim, msh_side.nqn, 1, msh_side.nel]);
      normal = geopdes_prod__ (Jinv, normal);
      normal = reshape (normal, [msh.rdim, msh_side.nqn, msh_side.nel]);
      norms = reshape (geopdes_norm__ (normal), msh_side.nqn, msh_side.nel);
      msh_side.charlen = (2./norms);
      clear normal norms
    end
  end

  clear normal norms

end
