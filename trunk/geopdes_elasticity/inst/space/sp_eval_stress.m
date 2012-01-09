% SP_EVAL_STRESS: Compute the stress of the solution, given by its degrees of freedom, at a given set of points.
%
%   [eu, F] = sp_eval_stress (u, space, geometry, lambda, mu, pts);
%   [eu, F] = sp_eval_stress (u, space, geometry, lambda, mu, npts);
%
% INPUT:
%     
%     u:          vector of dof weights
%     space:      object defining the discrete space (see sp_bspline_2d)
%     geometry:   geometry structure (see geo_load)
%     lambda, mu: Lame' parameters
%     pts:        cell array with coordinates of points along each parametric direction
%     npts:       number of points along each parametric direction
%
% OUTPUT:
%
%     eu: the function evaluated at the given points 
%     F:  grid points in the physical domain, that is, the mapped points
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2012 Rafael Vazquez
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

function [stress, F] = sp_eval_stress (u, space, geometry, npts, lambda, mu)

  if (space.ncomp == 1)
    error ('sp_eval_stress: the stress is not computed for scalars')
  end

  ndim = numel (npts);

  if (iscell (npts))
    pts = npts;
    npts = cellfun (@numel, pts);
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
    msh = msh_2d ({brk{1}, brk{2}}, pts, [], geometry, 'boundary', false);
  elseif (ndim == 3)
    msh = msh_3d ({brk{1}, brk{2}, brk{3}}, pts, [], geometry, 'boundary', false);
  end
  sp  = space.constructor (msh);

  [stress, F] = sp_eval_stress_msh (u, sp, msh, lambda, mu);
  F  = reshape (F, [ndim, npts]);
  stress = reshape (stress, [sp.ncomp, ndim, npts]);

end
