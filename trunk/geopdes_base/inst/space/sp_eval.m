% SP_EVAL: Evaluate a function, given by its degrees of freedom, at a given set of points.
%
%   [eu, F] = sp_eval (u, space, geometry, pts);
%   [eu, F] = sp_eval (u, space, geometry, npts);
%
% INPUT:
%     
%     u:         vector of dof weights
%     space:     class defining the space (see sp_bspline_2d)
%     geometry:  geometry structure (see geo_load)
%     pts:       cell array with coordinates of points along each parametric direction
%     npts:      number of points along each parametric direction
%
% OUTPUT:
%
%     eu: the function evaluated in the given points 
%     F:  grid points in the physical domain, that is, the mapped points
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

function [eu, F] = sp_eval (u, space, geometry, npts);

  ndim = numel (npts);

  if (iscell (npts))
    pts = npts;
  elseif (isvector (npts))
    if (ndim == 2)
      pts = {(linspace (0, 1, npts(1))), (linspace (0, 1, npts(2)))};
    elseif (ndim == 3)
      pts = {(linspace (0, 1, npts(1))), (linspace (0, 1, npts(2))), (linspace (0, 1, npts(3)))};
    end
  end

  for jj = 1:ndim
    pts{jj} = pts{jj}(:)';
    if (numel (pts{jj}) > 1)
      brk{jj} = [0, pts{jj}(1:end-1) + diff(pts{jj})/2, 1]; 
    else
      brk{jj} = [0 1];
    end
  end

  if (ndim == 2)
    msh = msh_2d ({brk{1}, brk{2}}, pts, [], geometry, 'no boundary');
  elseif (ndim == 3)
    msh = msh_3d ({brk{1}, brk{2}, brk{3}}, pts, [], geometry, 'no boundary');
  end
  sp  = space.constructor (msh);

  [eu, F] = sp_eval_msh (u, sp, msh);

  if (ndim == 2)
    F  = reshape (F, ndim, numel (pts{1}), numel (pts{2}));
    eu = squeeze (reshape (eu, sp.ncomp, numel (pts{1}), numel (pts{2})));
  elseif (ndim == 3)
    F  = reshape (F, ndim, numel (pts{1}), numel (pts{2}), numel (pts{3}));
    eu = squeeze (reshape (eu, sp.ncomp, numel (pts{1}), numel (pts{2}), numel (pts{3})));
  end

end
