% MSH_EVAL_BOUNDARY_SIDE: evaluate the parameterization in one boundary side of the domain.
%
%     [msh_side, msh_side_from_interior] = msh_eval_boundary_side (msh, iside, [element_list]);
%
% INPUTS:
%     
%    msh:   mesh object (see msh_cartesian)
%    iside: number of the boundary side to compute, from 1 to 2*msh.ndim (see the file geo_specs for documentation about face numbering)
%    element_list: elements on which to compute the parametrization
%
% OUTPUT:
%
%     msh_side: structure that contains the following fields
%
%     FIELD_NAME    (SIZE)                  DESCRIPTION
%     side_number   (scalar)                  number of the side
%     nel           (scalar)                  number of elements of the boundary side
%     nel_dir       (1 x ndim vector)         number of elements in each parametric direction
%     nqn           (scalar)                  number of quadrature nodes per element
%     nqn_dir       (1 x ndim vector)         number of quadrature nodes per element in each parametric direction
%     quad_nodes    (ndim x nqn x nel vector) coordinates of the quadrature nodes in parametric domain
%     quad_weights  (nqn x nel vector)        weights associated to the quadrature nodes
%     geo_map       (rdim x nqn x nel vector) physical coordinates of the quadrature nodes
%     geo_map_jac   (rdim x ndim x nqn x nel) Jacobian matrix of the map evaluated at the quadrature nodes
%     geo_map_der2  (rdim x ndim x ndim x nqn x nel) Hessian matrix of the map evaluated at the quadrature nodes
%     geo_map_der3  (rdim x ndim x ndim x ndim x nqn x nel) Third derivatives tensor of the map evaluated at the quadrature nodes
%     geo_map_der4  (rdim x ndim x ndim x ndim x ndim x nqn x nel) Fourth derivatives tensor of the map evaluated at the quadrature nodes
%     jacdet        (nqn x nel)               element of length, area, volume (if rdim = ndim, determinant of the Jacobian)
%
%     msh_side_from_interior: mesh structure that contains quadrature
%       points on the boundary, but which computes information from the
%       volumetric parametrization, like the derivative in the normal direction
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2014, 2015, 2017 Rafael Vazquez
% Copyright (C) 2017 Luca Coradello
% Copyright (C) 2023 Pablo Antolin
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

function [msh_side, msh_side_from_interior] = msh_eval_boundary_side (msh, iside, element_list)

if (nargin == 2)
  msh_side = struct (msh_precompute (msh.boundary(iside)));
  msh_side.side_number = iside;

  if (msh.ndim > 1)

    %    ind  = [2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind  = [2 2 1 1] in 2D;
    %    ind2 = [1 1 2 2 3 3] in 3D,                  %ind2 = [1 1 2 2] in 2D
    ind2 = ceil (iside/2);
    ind = setdiff (1:msh.ndim, ind2);

    % Compute the normal vector. This requires the derivative also in the
    % normal direction (the boundary manifold is not enough).
    qn(ind) = msh_side.qn;
    if (mod (iside, 2) == 0)
      qn{ind2} = msh.breaks{ind2}(end);
    else
      qn{ind2} = msh.breaks{ind2}(1);
    end

    % Auxiliary vector sizes, to use with reshape and permute
    reorder = @(x) x(:)';
    psize = reorder ([msh_side.nqn_dir; msh_side.nel_dir]);
    vorder = [1:2:msh_side.ndim*2, 2:2:msh_side.ndim*2]; % [1 3 5 2 4 6], for ndim = 3

    jac = feval (msh.map_der, cellfun (reorder, qn, 'UniformOutput', false));
    jac = reshape (jac, [msh.rdim, msh.ndim, psize]);
    jac = permute (jac, [1 2 vorder+2]);

    geo_map_jac = reshape(jac, msh.rdim, msh.ndim, msh_side.nqn, msh_side.nel);
    JinvT = geopdes_invT__ (geo_map_jac);
    JinvT = reshape (JinvT, [msh.rdim, msh.ndim, msh_side.nqn, msh_side.nel]);

    normal = zeros (msh.ndim, msh_side.nqn, 1, msh_side.nel);
    normal(ind2,:) = (-1)^iside;
    normal = geopdes_prod__ (JinvT, normal);
    normal = reshape (normal, [msh.rdim, msh_side.nqn, msh_side.nel]);
    % Now normalize
    norms = reshape (geopdes_norm__ (normal), [1, msh_side.nqn, msh_side.nel]);
    msh_side.normal = bsxfun (@rdivide, normal, norms);


    % Compute the characteristic length in the normal direction
    % This has to be computed using information from the interior
    qn = msh.qn; qw = msh.qw; nel_dir = msh.nel_dir;
    if (mod (iside, 2) == 0)
      qn{ind2} = qn{ind2}(:,end);
      if (~isempty(msh.qw))
        qw = qw{ind2}(:,end);
      else
        qw = [];
      end
    else
      qn{ind2} = qn{ind2}(:,1);
      if (~isempty(msh.qw))
        qw = qw{ind2}(:,1);
      else
        qw = [];
      end
    end
    nel_dir(ind2) = 1;

    reorder = @(x) x(:)';
    psize = reorder ([msh.nqn_dir; nel_dir]);
    vorder = [1:2:msh.ndim*2, 2:2:msh.ndim*2]; % [1 3 5 2 4 6], for ndim = 3

    jac = feval (msh.map_der, cellfun (reorder, qn, 'UniformOutput', false));
    jac = reshape (jac, [msh.rdim, msh.ndim, psize]);
    jac = permute (jac, [1 2 vorder+2]);
    jac = reshape (jac, [msh.rdim, msh.ndim, msh.nqn_dir, msh_side.nel]);

    jac = reshape (sqrt (sum (jac(:,ind2,:,:,:,:).^2, 1)), [msh.nqn_dir, msh_side.nel]); % Module

    qsize = ones (1, msh.ndim);
    qsize(ind2) = numel (qw);
    jac_times_qw = bsxfun(@times, jac, reshape (qw, qsize));
    msh_side.charlen = reshape (sum (jac_times_qw, ind2), msh_side.nqn, msh_side.nel);
  end

elseif (nargin == 3)
    
  if (isempty (element_list))
    msh_side.quad_weights = [];
    msh_side.geo_map = [];
    msh_side.geo_map_jac = [];
    msh_side.geo_map_der2 = [];
    msh_side.geo_map_der3 = [];
    msh_side.geo_map_der4 = [];
    msh_side.jacdet = [];
    msh_side.element_size = [];
    msh_side.normal = [];
    return
  end

  msh_side = msh_evaluate_element_list (msh.boundary(iside), element_list);
  msh_side.side_number = iside;

  if (msh.ndim > 1)

    %    ind  = [2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind  = [2 2 1 1] in 2D;
    %    ind2 = [1 1 2 2 3 3] in 3D,                  %ind2 = [1 1 2 2] in 2D
    ind2 = ceil (iside/2);
    ind = setdiff (1:msh.ndim, ind2);

    % Compute the normal vector. This requires the derivative also in the
    % normal direction (the boundary manifold is not enough).
    element_list = element_list(:)';
    indices = cell (msh_side.ndim, 1);
    [indices{:}] = ind2sub (msh_side.nel_dir, element_list);
    indices = cell2mat (indices);
    
    qn_elems = arrayfun(@(ii) {msh.boundary(iside).qn{ii}(:,indices(ii,:))}, 1:msh_side.ndim);
    qqn = cell (1,msh_side.nel);
    for iel = 1:numel(element_list)
      for idim = 1:numel(ind)
        qqn{iel}{ind(idim)} = qn_elems{idim}(:,iel)';
      end
      if (mod (iside, 2) == 0)
        qqn{iel}{ind2} = msh.breaks{ind2}(end);
      else
        qqn{iel}{ind2} = msh.breaks{ind2}(1);
      end
    end    
    
    jac = cellfun (@(x) feval (msh.map_der, x), qqn, 'UniformOutput', false);
    geo_map_jac  = zeros (msh.rdim, msh.ndim, msh_side.nqn, msh_side.nel);
    for iel = 1:numel(element_list)
      geo_map_jac(:,:,:,iel) = jac{iel};
    end
    JinvT = geopdes_invT__ (geo_map_jac);
    JinvT = reshape (JinvT, [msh.rdim, msh.ndim, msh_side.nqn, msh_side.nel]);
    normal = zeros (msh.ndim, msh_side.nqn, 1, msh_side.nel);
    normal(ind2,:) = (-1)^iside;
    normal = geopdes_prod__ (JinvT, normal);
    normal = reshape (normal, [msh.rdim, msh_side.nqn, msh_side.nel]);
    % Now normalize
    norms = reshape (geopdes_norm__ (normal), [1, msh_side.nqn, msh_side.nel]);
    msh_side.normal = bsxfun (@rdivide, normal, norms);
   
  end
    
end

if (nargout == 2)
  brk_bnd = msh.breaks; qn_bnd = msh.qn; qw_bnd = msh.qw;
  ind2 = ceil (iside/2);
  if (mod (iside, 2) == 1)
    brk_bnd{ind2} = brk_bnd{ind2}(1:2);
    qn_bnd{ind2} = brk_bnd{ind2}(1);
    qw_bnd{ind2} = 1;
  else
    brk_bnd{ind2} = brk_bnd{ind2}(end-1:end);
    qn_bnd{ind2} = brk_bnd{ind2}(end);
    qw_bnd{ind2} = 1;
  end

  geo.map = msh.map; geo.map_der = msh.map_der; geo.map_der2 = msh.map_der2; geo.map_der3 = msh.map_der3; geo.map_der4 = msh.map_der4;
  geo.rdim = msh.rdim;
  msh_side_from_interior = msh_cartesian (brk_bnd, qn_bnd, qw_bnd, geo, 'boundary', false);
  if (nargin == 2)
    msh_side_from_interior = msh_precompute (msh_side_from_interior);
  elseif (nargin == 3)
    msh_side_from_interior = msh_evaluate_element_list (msh_side_from_interior, element_list);
  end
end


end
