% MSH_3D_TENSOR_PRODUCT: construct a 3d tensor product mesh
%
%     msh = msh_3d_tensor_product(breaks, qn, qw, [opts ...]);
%
% INPUTS:
%     
%     breaks:  breaks along each direction in parametric space (repetitions are ignored)
%     qn:      quadrature nodes along each direction in parametric space
%     qw:      quadrature weights along each direction in parametric space
%     block:   which elements are assigned to this processor
%     opts:    'no boundary', no boundary mesh is built
%              'boundary' followed a vector containing the indices of
%              the nedded boundaries. Boundaries are labelled this way
%                   ________ 
%                  /   4   /|
%                 / ______/ |
%               1 |       |2|
%                 |   5   | /
%                 |_______|/ 
%                     3
%
%              3 is the bottom face and 6 the back one.
% OUTPUT:
%
%     msh.nel:          number of elements
%     msh.breaks:       unique(breaks)
%     msh.nqn:          number of nodes per element
%     msh.quad_nodes:   node coordinates in parametric space
%     msh.quad_weights: quadrature weights
% 
% Copyright (C) 2009, 2010 Carlo de Falco
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 2 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function msh = msh_3d_tensor_product (breaks, qn, qw, varargin)

  needed_boundary=1:6; % by default all boundaries are computed

  if (~isempty (varargin))
    ii=1;
    while (ii<=length(varargin))
      if (strcmpi (varargin {ii}, 'no boundary'))
        needed_boundary=[];
      elseif (strcmpi (varargin {ii}, 'boundary'))
        if (ii+1 <= length(varargin))
          needed_boundary = sort( varargin {ii+1});
          ii=ii+1;
        else
          error ('msh_3d_tensor_product: ''boundary'' option must be passed in the [option, value] format');
        end
      else
        error ('msh_3d_tensor_product: unknown option %s', varargin {ii});
      end
      ii=ii+1;
    end
  end

  msh.qn = qn;

  % this only for precaution, breaks should already not have any repetitions
  msh.breaks = {unique(breaks{1}), unique(breaks{2}), unique(breaks{3})};

  nelu = numel (msh.breaks{1}) - 1;
  nelv = numel (msh.breaks{2}) - 1;
  nelw = numel (msh.breaks{3}) - 1;

  qnu = qn{1};  qnv = qn{2}; qnw = qn{3};
  nqnu = size (qnu, 1); nqnv = size (qnv, 1); nqnw = size (qnw, 1);

  nel  = nelu * nelv * nelw;
  msh.nel  = nel;
  msh.nqn  = nqnu * nqnv * nqnw;
  
  quad_nodes_u = reshape (qnu, nqnu, 1, 1, nelu, 1, 1);
  quad_nodes_u = repmat  (quad_nodes_u, [1, nqnv, nqnw, 1, nelv, nelw]);
  quad_nodes_u = reshape (quad_nodes_u, [], nel);

  quad_nodes_v = reshape (qnv, 1, nqnv, 1, 1, nelv, 1);
  quad_nodes_v = repmat  (quad_nodes_v, [nqnu, 1, nqnw, nelu, 1, nelw]);
  quad_nodes_v = reshape (quad_nodes_v, [], nel);

  quad_nodes_w = reshape (qnw, 1, 1, nqnw, 1, 1, nelw);
  quad_nodes_w = repmat  (quad_nodes_w, [nqnu, nqnv, 1, nelu, nelv, 1]);
  quad_nodes_w = reshape (quad_nodes_w, [], nel);

  msh.quad_nodes(1, :, :) = quad_nodes_u;
  msh.quad_nodes(2, :, :) = quad_nodes_v;
  msh.quad_nodes(3, :, :) = quad_nodes_w;

  clear quad_nodes_u quad_nodes_v quad_nodes_w

  if (~isempty (qw))
    qwu = qw{1};  qwv = qw{2}; qww = qw{3};
    quad_weights_u = reshape (qwu, nqnu, 1, 1, nelu, 1, 1);
    quad_weights_u = repmat  (quad_weights_u, [1, nqnv, nqnw, 1, nelv, nelw]);
    quad_weights_u = reshape (quad_weights_u, [], nel);

    quad_weights_v = reshape (qwv, 1, nqnv, 1, 1, nelv, 1);
    quad_weights_v = repmat  (quad_weights_v, [nqnu, 1, nqnw, nelu, 1, nelw]);
    quad_weights_v = reshape (quad_weights_v, [], nel);

    quad_weights_w = reshape (qww, 1, 1, nqnw, 1, 1, nelw);
    quad_weights_w = repmat  (quad_weights_w, [nqnu, nqnv, 1, nelu, nelv, 1]);
    quad_weights_w = reshape (quad_weights_w, [], nel);

    msh.quad_weights = quad_weights_u .* quad_weights_v .* quad_weights_w;

    edge = zeros(1,3);
    for idir=1:3
      edge(idir)=msh.breaks{idir}(end)-msh.breaks{idir}(1);
    end
    area = prod(edge);
    if (abs (sum ( msh.quad_weights(:)) - area) > 1e-10)
      warning ('msh_3d_tensor_product: inconsistent quadrature formula')
    end

    clear quad_weights_u quad_weights_v quad_weights_w
  end

  msh.geo_map = msh.quad_nodes;
  msh.geo_map_jac = zeros (3, 3, msh.nqn, msh.nel);
  msh.geo_map_jac(1, 1, :, :) = 1;
  msh.geo_map_jac(2, 2, :, :) = 1;
  msh.geo_map_jac(3, 3, :, :) = 1;
  msh.jacdet = ones (msh.nqn, msh.nel);
  
  msh.boundary_list = needed_boundary; 
  for iside = msh.boundary_list
%%  ind =[2 3; 2 3; 1 3; 1 3; 1 2; 1 2]
    ind = setdiff (1:3, ceil(iside/2)); 
%%  ind2 = [1 1 2 2 3 3];
    ind2 = floor ((iside+1)/2);
      
    bnd_aux = msh_2d_tensor_product (msh.breaks(ind), msh.qn(ind), qw(ind), 'no boundary');

    boundary = bnd_aux;

    if (mod (iside, 2) == 1)
      boundary.quad_nodes = ones ([3, bnd_aux.nqn bnd_aux.nel])*breaks{ind2}(1);   %% this can be put to 0
    else
      boundary.quad_nodes = ones ([3, bnd_aux.nqn bnd_aux.nel])*breaks{ind2}(end); %% this can be put to 1
    end
    
    boundary.quad_nodes(ind,:,:) = bnd_aux.quad_nodes;
    boundary.geo_map = boundary.quad_nodes;
    boundary.geo_map_jac = zeros (3, 3, boundary.nqn, boundary.nel);
    boundary.geo_map_jac (1, 1, : , :) = 1;
    boundary.geo_map_jac (2, 2, : , :) = 1;
    boundary.geo_map_jac (3, 3, : , :) = 1;
    boundary.jacdet = ones(boundary.nqn, boundary.nel);
    boundary.normal = zeros (3, boundary.nqn, boundary.nel);
    boundary.normal(ind2,:,:) = (-1)^iside;
      
    msh.boundary(iside) = boundary;
  end

end
