% MSH_3D: constructor of the class for 3d tensor product meshes.
%
%     msh = msh_3d (breaks, qn, qw, geometry, opts);
%
% INPUTS:
%     
%     breaks:   breaks along each direction in parametric space (repetitions are ignored)
%     qn:       quadrature nodes along each direction in parametric space
%     qw:       quadrature weights along each direction in parametric space
%     geometry: structure representing the geometrical mapping
%     opts:     if opts == 'no boundary', the boundary terms are not computed
%   
% OUTPUT:
%
%     msh: class containing the following fields and methods
%
%     FIELD_NAME    (SIZE)                  DESCRIPTION
%     nel           (scalar)                number of elements of the partition
%     nelu          (scalar)                number of elements in the first parametric direction
%     nelv          (scalar)                number of elements in the second parametric direction
%     nelw          (scalar)                number of elements in the third parametric direction
%     nelcol        (scalar)                number of elements in one "column" of the mesh (actually, nelv*nelw)
%     nqn           (scalar)                number of quadrature nodes per element
%     nqnu          (scalar)                number of quadrature nodes per element in the first parametric direction
%     nqnv          (scalar)                number of quadrature nodes per element in the second parametric direction
%     nqnw          (scalar)                number of quadrature nodes per element in the third parametric direction
%     breaks        (1 x 3 cell-array)      unique(breaks)
%     qn            (1 x 3 cell-array)      quadrature nodes along each direction in parametric domain
%     qw            (1 x 3 cell-array)      quadrature weights along each direction in parametric space
%     boundary      (1 x 6 struct-array)    it contains a one-dimensional 'msh' structure for each face of the boundary 
%
%     METHOD NAME
%     msh_evaluate_col: computes the parameterization (and its derivatives) of
%                       the quadrature points in one column of the mesh, i.e.,
%                       fixing the element in the first parametric direction.
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

function msh = msh_3d (breaks, qn, qw, geo, opts)

  if (nargin < 5) 
    opts = '';
  end

  msh.qn = qn;
  msh.qw = qw;

  % this only for precaution, breaks should already not have any repetitions
  msh.breaks = {unique(breaks{1}), unique(breaks{2}), unique(breaks{3})};

  msh.nelu = numel (msh.breaks{1}) - 1;
  msh.nelv = numel (msh.breaks{2}) - 1;
  msh.nelw = numel (msh.breaks{3}) - 1;
  msh.nel = msh.nelu * msh.nelv * msh.nelw;
  msh.nelcol = msh.nelv * msh.nelw;

  qnu = qn{1};  qnv = qn{2}; qnw = qn{3};
  msh.nqnu = size (qnu,1); 
  msh.nqnv = size (qnv,1);
  msh.nqnw = size (qnw,1);
  msh.nqn  = msh.nqnu * msh.nqnv * msh.nqnw;
  
  if (~strcmpi (opts, 'no boundary'))
    warning ('Ricordati di cambiare il bordo 3d')
    for iside = 1:6
%%    ind =[2 3; 2 3; 1 3; 1 3; 1 2; 1 2]
      ind = setdiff (1:3, ceil(iside/2)); 
%%    ind2 = [1 1 2 2 3 3];
      ind2 = floor ((iside+1)/2);

      bnd_aux = msh_2d_tensor_product (msh.breaks(ind), msh.qn(ind), qw(ind), 'no boundary');

      boundary = bnd_aux;

      if (mod (iside, 2) == 1)
        boundary.quad_nodes = zeros ([3, bnd_aux.nqn bnd_aux.nel]);
      else
        boundary.quad_nodes = ones ([3, bnd_aux.nqn bnd_aux.nel]);
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

      qnu = boundary.quad_nodes(1,:,:);
      qnv = boundary.quad_nodes(2,:,:);
      qnw = boundary.quad_nodes(3,:,:);
      F   = feval (geo.map, [qnu(:), qnv(:), qnw(:)]');
      jac = feval (geo.map_der, [qnu(:), qnv(:), qnw(:)]');

      boundary.geo_map = ...
           reshape (F, size (boundary.quad_nodes));
      boundary.geo_map_jac = ...
           reshape(jac, 3, 3, boundary.nqn, boundary.nel);

      jacdet = ...
        geopdes_norm__ (cross (squeeze (boundary.geo_map_jac(:,ind(1),:,:)), ...
                  squeeze (boundary.geo_map_jac(:,ind(2),:,:))));
      boundary.jacdet = reshape (jacdet, boundary.nqn, boundary.nel);

      [JinvT, jacdet] = geopdes_invT__ (boundary.geo_map_jac);
      JinvT = reshape (JinvT, [3, 3, boundary.nqn, boundary.nel]);

      normal = reshape (boundary.normal, [3, boundary.nqn, 1, boundary.nel]);
      normal = geopdes_prod__ (JinvT, normal);
      normal = reshape (normal, [3, boundary.nqn, boundary.nel]);
      norms = repmat (reshape (geopdes_norm__ (normal), [1, boundary.nqn, boundary.nel]), [3 1 1]);
      boundary.normal = normal ./ norms;

      msh.boundary(iside) = boundary;
    end
  else
    msh.boundary = [];    
  end

  msh.map = geo.map;
  msh.map_der = geo.map_der;
  msh = class (msh, 'msh_3d');

end
