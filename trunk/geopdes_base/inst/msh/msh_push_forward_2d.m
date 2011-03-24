% MSH_PUSH_FORWARD_2D: construct a msh structure by applying the geometry map to a msh structure in the reference domain.
%
%     msh = msh_push_forward_2d (msh, geo, 'option1', value1, ...)
%
% INPUTS:
%
%     msh:  parametric domain partition and quadrature nodes
%     geo:  structure representing the geometrical mapping
%     'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+-----------
%               der2    |      false      |  compute second order derivatives
%                       |                 |  of the geometry at quad nodes
%
% OUTPUT:
%
%     msh: structure containing the domain partition and the quadrature rule in both the parametric and the physical domain, which contains the following fields
%
%          FIELD_NAME    (SIZE)                  DESCRIPTION
%          nel           (scalar)                number of elements of the partition
%          nqn           (scalar)                number of quadrature nodes per element
%          breaks        (1 x 2 cell-array)      unique(breaks)
%          quad_nodes    (2 x nqn x nel vector)  coordinates of the quadrature nodes in parametric space
%          quad_weights  (nqn x nel vector)      weights associated to the quadrature nodes
%          geo_map       (2 x nqn x nel vector)  physical coordinates of the quadrature nodes
%          geo_map_jac   (2 x 2 x nqn x nel)     Jacobian matrix of the map evaluated at the quadrature nodes
%          geo_map_der2  (2 x 2 x 2 x nqn x nel) Second order derivatives of the map evaluated at the quadrature nodes ()
%          jacdet        (nqn x nel)             determinant of the Jacobian evaluated in the quadrature points
%          boundary      (1 x 4 struct-array)    it contains a one-dimensional 'msh' structure for each edge of the boundary 
%  For more details, see the documentation
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

function msh = msh_push_forward_2d (msh, geo, varargin)

der2 = false;
if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('msh_push_forward_2d: options must be passed in the [option, value] format');
  end
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'der2'))
      der2 = varargin {ii+1};
    else
      error ('msh_push_forward_2d: unknown option %s', varargin {ii});
    end
  end
end

qnu = msh.quad_nodes(1,:,:);
qnv = msh.quad_nodes(2,:,:);
F   = feval (geo.map, [qnu(:), qnv(:)]');
jac = feval (geo.map_der, [qnu(:), qnv(:)]');

if (der2)
  if (isfield (geo, 'map_der2'))
    msh.geo_map_der2 = reshape (feval (geo.map_der2, [qnu(:), qnv(:)]'), 2, 2, 2, msh.nqn, msh.nel);
  else 
    error ('msh_push_forward_2d: a function to compute second order derivatives has not been provided')
  end
end

msh.geo_map = reshape (F, size (msh.quad_nodes));
msh.geo_map_jac = reshape (jac, 2, 2, msh.nqn, msh.nel);
msh.jacdet = abs (geopdes_det__ (msh.geo_map_jac));
msh.jacdet = reshape (msh.jacdet, [msh.nqn, msh.nel]);


if (isfield (msh, 'boundary'))
  for iside = 1:4
    ind1 = floor ((iside+1)/2);  %ind1 = [1 1 2 2];
    ind2 = mod (ind1, 2) + 1;    %ind2 = [2 2 1 1];
    qnu = msh.boundary(iside).quad_nodes(1,:,:);
    qnv = msh.boundary(iside).quad_nodes(2,:,:);
    F   = feval (geo.map, [qnu(:), qnv(:)]');
    jac = feval (geo.map_der, [qnu(:), qnv(:)]');

    msh.boundary(iside).geo_map = ...
        reshape (F, size (msh.boundary(iside).quad_nodes));
    msh.boundary(iside).geo_map_jac = ...
        reshape(jac, 2, 2, msh.boundary(iside).nqn, msh.boundary(iside).nel);
    jacdet = ...
        geopdes_norm__ (squeeze (msh.boundary(iside).geo_map_jac(:,ind2,:,:)));
    msh.boundary(iside).jacdet = reshape (jacdet, ...
                                          msh.boundary(iside).nqn, msh.boundary(iside).nel);

    [JinvT, jacdet] = geopdes_invT__ (msh.boundary(iside).geo_map_jac);
    JinvT = reshape (JinvT, [2, 2, msh.boundary(iside).nqn, msh.boundary(iside).nel]);
    normal = reshape (msh.boundary(iside).normal, [2, msh.boundary(iside).nqn, 1, msh.boundary(iside).nel]);

    normal = geopdes_prod__ (JinvT, normal);
    normal = reshape (normal, [2, msh.boundary(iside).nqn, msh.boundary(iside).nel]);
    norms = repmat (reshape (geopdes_norm__ (normal), [1, msh.boundary(iside).nqn, msh.boundary(iside).nel]), [2 1 1]);
    msh.boundary(iside).normal = normal ./ norms;

  end
end

end
