% GEO_3D_NURBS: construct a geometry map from a structure of the NURBS toolbox.
%
%   output = geo_3d_nurbs (nurbs, pts, ders)
%
% INPUTS:
%
%   nurbs:  NURBS structure that defines the geometry
%   pts  :  points where the map has to be evaluated
%   ders :  number of derivatives to be evaluated (from 0 to 2)
%
% OUTPUT:
%
%   output: for ders = 0, the parametrization F evaluated at pts
%           for ders = 1, the Jacobian of the parametrization, evaluated at pts
%
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

function varargout = geo_3d_nurbs (nurbs, pts, ders)


  switch (ders)
    case 0
      F = nrbeval (nurbs, pts);
      varargout{1} = F;
    case 1
      deriv = nrbderiv (nurbs);
      [F, jac] = nrbdeval (nurbs, deriv, pts);
      if (iscell (pts))
        npts = prod (cellfun (@numel, pts));
        map_jac = zeros (3, 3, npts);
        map_jac(1:3, 1, :) = reshape (jac{1}, 3, 1, npts);
        map_jac(1:3, 2, :) = reshape (jac{2}, 3, 1, npts);
        map_jac(1:3, 3, :) = reshape (jac{3}, 3, 1, npts);
      else
        map_jac = zeros (3, 3, size (pts, 2));
        map_jac(1:3, 1, :) = reshape (jac{1}(1:3, :), 3, 1, size (pts, 2));
        map_jac(1:3, 2, :) = reshape (jac{2}(1:3, :), 3, 1, size (pts, 2));
        map_jac(1:3, 3, :) = reshape (jac{3}(1:3, :), 3, 1, size (pts, 2));
      end
      varargout{1} = map_jac;
    otherwise
      error ('geo_3d_nurbs: number of derivatives limited to 1 (until now)')
  end

end
