% MSH_2D_TENSOR_PRODUCT: construct a 2d tensor product mesh in the parametric domain.
%
%     msh = msh_2d_tensor_product (breaks, qn, qw, opts);
%
% INPUTS:
%     
%     breaks:  breaks along each direction in parametric space (repetitions are ignored)
%     qn:      quadrature nodes along each direction in parametric space
%     qw:      quadrature weights along each direction in parametric space
%     opts:    if opts == 'no boundary', the boundary terms are not computed
%   
% OUTPUT:
%
%     msh: structure containing the following fields
%
%          FIELD_NAME    (SIZE)                  DESCRIPTION
%          nel           (scalar)                number of elements of the partition
%          nqn           (scalar)                number of quadrature nodes per element
%          breaks        (1 x 2 cell-array)      unique(breaks)
%          quad_nodes    (2 x nqn x nel vector)  coordinates of the quadrature nodes in parametric space
%          quad_weights  (nqn x nel vector)      weights associated to the quadrature nodes
%          geo_map       (2 x nqn x nel vector)  physical coordinates of the quadrature nodes
%          geo_map_jac   (2 x 2 x nqn x nel)     Jacobian matrix of the map evaluated at the quadrature nodes
%          jacdet        (nqn x nel)             determinant of the Jacobian evaluated in the quadrature points
%          boundary      (1 x 4 struct-array)    it contains a one-dimensional 'msh' structure for each edge of the boundary 
%
% For more details, see the documentation
% 
% Copyright (C) 2009, 2010 Carlo de Falco
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

function msh = msh_2d_tensor_product (breaks, qn, qw, opts)

  if (nargin < 4) 
    opts = '';
  end

  msh.qn = qn;

  msh.breaks = {unique(breaks{1}), unique(breaks{2})}; % this only for precaution, breaks
                                                       % should already not have any repetitions

  nelu = numel (msh.breaks{1}) - 1;
  nelv = numel (msh.breaks{2}) - 1;

  qnu = qn{1};  qnv = qn{2};
  nqnu = size (qnu,1); nqnv = size (qnv,1);

  nel  = nelu * nelv;
  msh.nel  = nel;
  msh.nqn  = nqnu * nqnv;
  
  quad_nodes_u = reshape (qnu, nqnu, 1, nelu, 1);
  quad_nodes_u = repmat  (quad_nodes_u, [1, nqnv, 1, nelv]);
  quad_nodes_u = reshape (quad_nodes_u, [], nel);

  quad_nodes_v = reshape (qnv, 1, nqnv, 1, nelv);
  quad_nodes_v = repmat  (quad_nodes_v, [nqnu, 1, nelu, 1]);
  quad_nodes_v = reshape (quad_nodes_v, [], nel);

  msh.quad_nodes(1, :, :) = quad_nodes_u;
  msh.quad_nodes(2, :, :) = quad_nodes_v;

  clear quad_nodes_u quad_nodes_v

  if (~isempty (qw))
    qwu = qw{1};  qwv = qw{2};
    quad_weights_u = reshape (qwu, nqnu, 1, nelu, 1);
    quad_weights_u = repmat  (quad_weights_u, [1, nqnv, 1, nelv]);
    quad_weights_u = reshape (quad_weights_u, [], nel);

    quad_weights_v = reshape (qwv, 1, nqnv, 1, nelv);
    quad_weights_v = repmat  (quad_weights_v, [nqnu, 1, nelu, 1]);
    quad_weights_v = reshape (quad_weights_v, [], nel);

    msh.quad_weights = quad_weights_u .* quad_weights_v;
    if (abs (sum ( msh.quad_weights(:)) - 1) > 1e-10)
      warning ('msh_2d_tensor_product: inconsistent quadrature formula')
    end

    clear quad_weights_u quad_weights_v
  end

  msh.geo_map = msh.quad_nodes;
  msh.geo_map_jac = zeros (2, 2, msh.nqn, msh.nel);
  msh.geo_map_jac(1, 1, :, :) = 1;
  msh.geo_map_jac(2, 2, :, :) = 1;
  msh.jacdet = ones (msh.nqn, msh.nel);


  if (~strcmpi (opts, 'no boundary'))
    for iside = 1:4
      ind = mod (floor ((iside+1)/2), 2) + 1;  %ind = [2 2 1 1];
      ind2 = floor ((iside+1)/2);              %ind2 = [1 1 2 2];
      boundary.breaks = msh.breaks{ind};
      boundary.nel = numel (msh.breaks{ind}) - 1;
      boundary.nqn = size (qn{ind},1);
      boundary.quad_weights = qw{ind};
      boundary.quad_nodes = zeros ([2, size(qn{ind})]);
      boundary.quad_nodes(ind,:,:) = qn{ind};
      if (iside == 2 || iside == 4)
        boundary.quad_nodes(ind2,:,:) = 1;
      end

      boundary.geo_map = boundary.quad_nodes;
      boundary.geo_map_jac = zeros (2, 2, boundary.nqn, boundary.nel);
      boundary.geo_map_jac (1, 1, :, :) = 1;
      boundary.geo_map_jac (2, 2, :, :) = 1;
      boundary.jacdet = ones (boundary.nqn, boundary.nel);

      boundary.normal = zeros (2, boundary.nqn, boundary.nel);
      boundary.normal(ind2,:,:) = (-1)^iside;

      msh.boundary(iside) = boundary;
    end
  end

end
