% GEO_DEFORM: create the geometry structure for a deformed NURBS entity.
%
%   new_geometry = geo_deform (u, space, geometry);
%
% INPUT:
%     
%     u:         vector of dof weights for the displacement
%     space:     structure representing the space of discrete functions
%     geometry:  geometry structure. It must contain a nurbs substructure
%                 as given in the NURBS toolbox.
%
% OUTPUT:
%
%     new_geometry: geometry structure for the deformed domain.
% 
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

function new_geom = geo_deform (u, space, geometry);

  if (~isfield (geometry, 'nurbs'))
    error ('geo_deform: only NURBS-based geometries are allowed')
  end

  if (space.ncomp == 2)
    nurbs = geometry.nurbs;
    ndof_comp = prod (space.ndof_dir, 2);

    u1 = reshape (u(1:ndof_comp(1)), [1, space.ndof_dir(1,:)]);
    u2 = reshape (u(ndof_comp(1)+[1:ndof_comp(2)]), [1, space.ndof_dir(2,:)]);
    weights = nurbs.coefs(4,:,:);

    nurbs.coefs(1,:,:) = nurbs.coefs(1,:,:) + u1 .* weights;
    nurbs.coefs(2,:,:) = nurbs.coefs(2,:,:) + u2 .* weights;

    new_geom.nurbs = nurbs;

    new_geom.map      = @(PTS) geo_2d_nurbs (new_geom.nurbs, PTS, 0);
    new_geom.map_der  = @(PTS) geo_2d_nurbs (new_geom.nurbs, PTS, 1);
    new_geom.map_der2 = @(PTS) geo_2d_nurbs (new_geom.nurbs, PTS, 2);

  elseif (space.ncomp == 3)
    nurbs = geometry.nurbs;
    ndof_comp = prod (space.ndof_dir, 2);

    u1 = reshape (u(1:ndof_comp(1)), [1, space.ndof_dir(1,:)]);
    u2 = reshape (u(ndof_comp(1)+[1:ndof_comp(2)]), [1, space.ndof_dir(2,:)]);
    u3 = reshape (u(ndof_comp(1)+ndof_comp(2)+[1:ndof_comp(3)]), [1, space.ndof_dir(3,:)]);
    weights = nurbs.coefs(4,:,:,:);

    nurbs.coefs(1,:,:,:) = nurbs.coefs(1,:,:,:) + u1 .* weights;
    nurbs.coefs(2,:,:,:) = nurbs.coefs(2,:,:,:) + u2 .* weights;
    nurbs.coefs(3,:,:,:) = nurbs.coefs(3,:,:,:) + u3 .* weights;

    new_geom.nurbs = nurbs;

    new_geom.map     = @(PTS) geo_3d_nurbs (new_geom.nurbs, PTS, 0);
    new_geom.map_der = @(PTS) geo_3d_nurbs (new_geom.nurbs, PTS, 1);

  end

end
