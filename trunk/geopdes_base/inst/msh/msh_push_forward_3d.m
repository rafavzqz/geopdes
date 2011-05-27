% MSH_PUSH_FORWARD_3D: construct a msh structure by applying the geometry map to a msh structure in the reference domain.
%
%     msh = msh_push_forward_3d (msh, geo)
%
% INPUTS:
%     
%     msh:  parametric domain partition and quadrature nodes
%   
% OUTPUT:
%
%     msh: structure containing the domain partition and the quadrature rule in both the parametric and the physical domain, which contains the following fields
%          FIELD_NAME    (SIZE)                  DESCRIPTION
%          nel           (scalar)                number of elements of the partition
%          nqn           (scalar)                number of quadrature nodes per element
%          breaks        (1 x 3 cell-array)      unique(breaks)
%          quad_nodes    (3 x nqn x nel vector)  coordinates of the quadrature nodes in parametric space
%          quad_weights  (nqn x nel vector)      weights associated to the quadrature nodes
%          geo_map       (3 x nqn x nel vector)  physical coordinates of the quadrature nodes
%          geo_map_jac   (3 x 3 x nqn x nel)     Jacobian matrix of the map evaluated at the quadrature nodes
%          jacdet        (nqn x nel)             determinant of the Jacobian evaluated at the quadrature points
%          boundary      (1 x 6 struct-array)    it contains a two-dimensional 'msh' structure for each face of the boundary 
%
%  For more details, see the documentation
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

function msh = msh_push_forward_3d (msh, geo)

  if (isfield (msh, 'qn'))
    qn = msh.qn;
    F   = feval (geo.map, {qn{1}(:)', qn{2}(:)' qn{3}(:)'});
    F = reshape (F, [3, msh.nqnu, msh.nelu, msh.nqnv, msh.nelv, msh.nqnw, msh.nelw]);
    F = permute (F, [1 2 4 6 3 5 7]);
    msh.geo_map = reshape (F, [3, msh.nqn, msh.nel]);

    jac = feval (geo.map_der, {qn{1}(:)', qn{2}(:)' qn{3}(:)'});
    jac = reshape (jac, [3, 3, msh.nqnu, msh.nelu, msh.nqnv, msh.nelv, msh.nqnw, msh.nelw]);
    jac = permute (jac, [1 2 3 5 7 4 6 8]);
    msh.geo_map_jac = reshape (jac, 3, 3, msh.nqn, msh.nel);
    msh.jacdet = abs (geopdes_det__ (msh.geo_map_jac));
    msh.jacdet = reshape (msh.jacdet, [msh.nqn, msh.nel]);
  else
    qnu = msh.quad_nodes(1,:,:);
    qnv = msh.quad_nodes(2,:,:);
    qnw = msh.quad_nodes(3,:,:);
    F   = feval (geo.map, [qnu(:), qnv(:) qnw(:)]');
    jac = feval (geo.map_der, [qnu(:), qnv(:) qnw(:)]');

    msh.geo_map = reshape (F, size (msh.quad_nodes));
    msh.geo_map_jac = reshape (jac, 3, 3, msh.nqn, msh.nel);
    msh.jacdet = abs (geopdes_det__ (msh.geo_map_jac));
    msh.jacdet = reshape (msh.jacdet, [msh.nqn, msh.nel]);
 end

  if (isfield (msh, 'boundary'))
    for iside = 1:6
      qnu = msh.boundary(iside).quad_nodes(1,:,:);
      qnv = msh.boundary(iside).quad_nodes(2,:,:);
      qnw = msh.boundary(iside).quad_nodes(3,:,:);
      F   = feval (geo.map, [qnu(:), qnv(:), qnw(:)]');
      jac = feval (geo.map_der, [qnu(:), qnv(:), qnw(:)]');

      msh.boundary(iside).geo_map = ...
           reshape (F, size (msh.boundary(iside).quad_nodes));
      msh.boundary(iside).geo_map_jac = ...
           reshape(jac, 3, 3, msh.boundary(iside).nqn, msh.boundary(iside).nel);

      switch (iside)
        case {1, 2}
          ind1 = 1;
          ind2 = [2, 3];
        case {3, 4}
          ind1 = 2;
          ind2 = [1, 3];
        case {5, 6}
          ind1 = 3;
          ind2 = [1, 2];
      end
      jacdet = ...
        geopdes_norm__ (cross (squeeze (msh.boundary(iside).geo_map_jac(:,ind2(1),:,:)), ...
                  squeeze (msh.boundary(iside).geo_map_jac(:,ind2(2),:,:))));
      msh.boundary(iside).jacdet = reshape (jacdet, ...
           msh.boundary(iside).nqn, msh.boundary(iside).nel);

      [JinvT, jacdet] = geopdes_invT__ (msh.boundary(iside).geo_map_jac);
      JinvT = reshape (JinvT, [3, 3, msh.boundary(iside).nqn, msh.boundary(iside).nel]);

      normal = reshape (msh.boundary(iside).normal, [3, msh.boundary(iside).nqn, 1, msh.boundary(iside).nel]);
      normal = geopdes_prod__ (JinvT, normal);
      normal = reshape (normal, [3, msh.boundary(iside).nqn, msh.boundary(iside).nel]);
      norms = repmat (reshape (geopdes_norm__ (normal), [1, msh.boundary(iside).nqn, msh.boundary(iside).nel]), [3 1 1]);
      msh.boundary(iside).normal = normal ./ norms;

    end
  end

end
