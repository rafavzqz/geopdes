% MSH_2D: constructor of the class for 2d tensor product meshes.
%
%     msh = msh_2d (breaks, qn, qw, geometry, opts);
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
%     nelcol        (scalar)                number of elements in one "column" of the mesh (actually, nelv)
%     nqn           (scalar)                number of quadrature nodes per element
%     nqnu          (scalar)                number of quadrature nodes per element in the first parametric direction
%     nqnv          (scalar)                number of quadrature nodes per element in the second parametric direction
%     breaks        (1 x 2 cell-array)      unique(breaks)
%     qn            (1 x 2 cell-array)      quadrature nodes along each direction in parametric domain
%     qw            (1 x 2 cell-array)      quadrature weights along each direction in parametric space
%     boundary      (1 x 4 struct-array)    it contains a one-dimensional 'msh' structure for each edge of the boundary 
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

function msh = msh_2d (breaks, qn, qw, geo, opts)

  if (nargin < 5) 
    opts = '';
  end

  msh.qn = qn;
  msh.qw = qw;

  msh.breaks = {unique(breaks{1}), unique(breaks{2})}; % this only for precaution, breaks
                                                       % should already not have any repetitions

  msh.nelu = numel (msh.breaks{1}) - 1;
  msh.nelv = numel (msh.breaks{2}) - 1;
  msh.nel = msh.nelu * msh.nelv;
  msh.nelcol = msh.nelv;

  qnu = qn{1};  qnv = qn{2};
  msh.nqnu = size (qnu,1); 
  msh.nqnv = size (qnv,1);
  msh.nqn  = msh.nqnu * msh.nqnv;
  
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

      qn1 = boundary.quad_nodes(1,:,:);
      qn2 = boundary.quad_nodes(2,:,:);
      F   = feval (geo.map, [qn1(:), qn2(:)]');
      jac = feval (geo.map_der, [qn1(:), qn2(:)]');

      boundary.geo_map = reshape (F, size (boundary.quad_nodes));
      boundary.geo_map_jac = reshape(jac, 2, 2, boundary.nqn, boundary.nel);
      boundary.jacdet = geopdes_norm__ (squeeze (boundary.geo_map_jac(:,ind,:,:)));

      normal = zeros (2, boundary.nqn, boundary.nel);
      normal(ind2,:,:) = (-1)^iside;

      [JinvT, jacdet] = geopdes_invT__ (boundary.geo_map_jac);
      JinvT = reshape (JinvT, [2, 2, boundary.nqn, boundary.nel]);

      normal = reshape (normal, [2, boundary.nqn, 1, boundary.nel]);
      normal = geopdes_prod__ (JinvT, normal);
      normal = reshape (normal, [2, boundary.nqn, boundary.nel]);
      norms = repmat (reshape (geopdes_norm__ (normal), [1, boundary.nqn, boundary.nel]), [2 1 1]);
      boundary.normal = normal ./ norms;

      msh.boundary(iside) = boundary;
    end
  else
    msh.boundary = [];    
  end

  msh.map = geo.map;
  msh.map_der = geo.map_der;
  msh = class (msh, 'msh_2d');

end
