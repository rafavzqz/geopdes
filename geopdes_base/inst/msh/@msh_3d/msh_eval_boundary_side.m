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
%     jacdet        (nqn x nel)             determinant of the Jacobian evaluated at the quadrature points
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

function msh_side = msh_eval_boundary_side (msh, iside)

  ind = setdiff (1:3, ceil(iside/2)); % ind = [2 3; 2 3; 1 3; 1 3; 1 2; 1 2] 
  ind2 = floor ((iside+1)/2); % ind2 = [1 1 2 2 3 3];

  msh_side = msh.boundary(iside);

  quad_weights_u = reshape (msh_side.qw{1}, msh_side.nqn_dir(1), 1, msh_side.nel_dir(1), 1);
  quad_weights_u = repmat  (quad_weights_u, [1, msh_side.nqn_dir(2), 1, msh_side.nel_dir(2)]);
  quad_weights_u = reshape (quad_weights_u, [], msh_side.nel);

  quad_weights_v = reshape (msh_side.qw{2}, 1, msh_side.nqn_dir(2), 1, msh_side.nel_dir(2));
  quad_weights_v = repmat  (quad_weights_v, [msh_side.nqn_dir(1), 1, msh_side.nel_dir(1), 1]);
  quad_weights_v = reshape (quad_weights_v, [], msh_side.nel);

  msh_side.quad_weights = quad_weights_u .* quad_weights_v;
  clear quad_weights_u quad_weights_v

  quad_nodes_u = reshape (msh_side.qn{1}, msh_side.nqn_dir(1), 1, msh_side.nel_dir(1), 1);
  quad_nodes_u = repmat  (quad_nodes_u, [1, msh_side.nqn_dir(2), 1, msh_side.nel_dir(2)]);
  quad_nodes_u = reshape (quad_nodes_u, [], msh_side.nel);

  quad_nodes_v = reshape (msh_side.qn{2}, 1, msh_side.nqn_dir(2), 1, msh_side.nel_dir(2));
  quad_nodes_v = repmat  (quad_nodes_v, [msh_side.nqn_dir(1), 1, msh_side.nel_dir(1), 1]);
  quad_nodes_v = reshape (quad_nodes_v, [], msh_side.nel);

  quad_nodes(1, :, :) = quad_nodes_u;
  quad_nodes(2, :, :) = quad_nodes_v;
  clear quad_nodes_u quad_nodes_v

  if (mod (iside, 2) == 0)
    msh_side.quad_nodes = ones ([3, msh_side.nqn, msh_side.nel]);
  else
    msh_side.quad_nodes = zeros ([3, msh_side.nqn, msh_side.nel]);
  end
  msh_side.quad_nodes(ind, :, :) = quad_nodes;

  qn1 = msh_side.quad_nodes(1,:,:);
  qn2 = msh_side.quad_nodes(2,:,:);
  qn3 = msh_side.quad_nodes(3,:,:);
  F   = feval (msh.map, [qn1(:), qn2(:), qn3(:)]');
  jac = feval (msh.map_der, [qn1(:), qn2(:), qn3(:)]');

  msh_side.geo_map = reshape (F, size (msh_side.quad_nodes));
  msh_side.geo_map_jac = reshape(jac, 3, 3, msh_side.nqn, msh_side.nel);

  jacdet = geopdes_norm__ (cross ...
                          (squeeze (msh_side.geo_map_jac(:,ind(1),:,:)), ...
                           squeeze (msh_side.geo_map_jac(:,ind(2),:,:))));
  msh_side.jacdet = reshape (jacdet, msh_side.nqn, msh_side.nel);

  JinvT = geopdes_invT__ (msh_side.geo_map_jac);
  JinvT = reshape (JinvT, [3, 3, msh_side.nqn, msh_side.nel]);

  normal = zeros (3, msh_side.nqn, 1, msh_side.nel);
  normal(ind2,:,:,:) = (-1)^iside;

  normal = geopdes_prod__ (JinvT, normal);
  normal = reshape (normal, [3, msh_side.nqn, msh_side.nel]);
  norms = repmat (reshape (geopdes_norm__ (normal), [1, msh_side.nqn, msh_side.nel]), [3 1 1]);
  msh_side.normal = normal ./ norms;
  clear normal norms

end
