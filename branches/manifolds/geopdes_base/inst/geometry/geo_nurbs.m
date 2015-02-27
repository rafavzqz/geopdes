% GEO_NURBS: construct a geometry map from a structure of the NURBS toolbox.
%
%   output = geo_nurbs (nurbs, pts, ders)
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
%           for ders = 2, the Hessian of the parametrization, evaluated at pts
%          
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011 Rafael Vazquez
% Copyright (C) 2013 Elena Bulgarello, Carlo de Falco, Sara Frizziero
% Copyright (C) 2015 Rafael Vazquez
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

function varargout = geo_nurbs (nurbs, pts, ders)

  ndim = numel (nurbs.order);

  if (any (abs(nurbs.coefs(3,:)) > 1e-12))
    rdim = 3;
  elseif (any (abs(nurbs.coefs(2,:)) > 1e-12))
    rdim = 2;
  else
    rdim = 1;
  end

  if (rdim > ndim)
    error ('geo_nurbs: the dimensions of your geometry seem to be wrong')
  end
  
  switch (ders)
    case 0
      F = nrbeval (nurbs, pts);
      varargout{1} = F(1:rdim, :);
    case 1
      deriv = nrbderiv (nurbs);
      [~, jac] = nrbdeval (nurbs, deriv, pts);
      if (iscell (pts))
        npts = prod (cellfun (@numel, pts));
        map_jac = zeros (rdim, ndim, npts);
        for idim = 1:ndim
          map_jac(1:rdim, idim, :) = reshape (jac{idim}(1:rdim,:,:), rdim, 1, npts);
        end
      else
        map_jac = zeros (rdim, ndim, size (pts, 2));
        for idim = 1:ndim
          map_jac(1:rdim, idim, :) = reshape (jac{idim}(1:rdim, :), rdim, 1, size (pts, 2));
        end
      end
      varargout{1} = map_jac;

    case 2
      [deriv, deriv2] = nrbderiv (nurbs);
      [~, ~, hessian] = nrbdeval (nurbs, deriv, deriv2, pts);

      if (iscell (pts))
        npts = prod (cellfun (@numel, pts));
        hess = zeros (rdim, ndim, ndim, npts);
        for idim = 1:ndim
          for jdim = 1:ndim
            hess(1:rdim, idim, jdim, :) = reshape (hessian{idim,jdim}(1:rdim,:,:), rdim, 1, 1, npts);
          end
        end
      else
        hess = zeros (rdim, ndim, ndim, size (pts, 2));
        for idim = 1:ndim
          for jdim = 1:ndim
            hess(1:rdim, idim, jdim, :) = reshape (hessian{idim,jdim}(1:rdim,:), rdim, 1, 1, size (pts, 2));
          end
        end
      end
      varargout{1} = hess;
    otherwise
      error ('geo_nurbs: number of derivatives limited to two')
  end

end
