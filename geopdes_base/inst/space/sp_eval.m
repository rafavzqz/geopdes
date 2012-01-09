% SP_EVAL: Compute the value or the derivatives of a function, given by its degrees of freedom, at a given set of points.
%
%   [eu, F] = sp_eval (u, space, geometry, pts, [option]);
%   [eu, F] = sp_eval (u, space, geometry, npts, [option]);
%
% INPUT:
%     
%     u:         vector of dof weights
%     space:     object defining the discrete space (see sp_bspline_2d)
%     geometry:  geometry structure (see geo_load)
%     pts:       cell array with coordinates of points along each parametric direction
%     npts:      number of points along each parametric direction
%     option:    accepted options are 'value' (default), 'gradient',
%                 and for vectors also 'curl', 'divergence'
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

function [eu, F] = sp_eval (u, space, geometry, npts, varargin)

  if (nargin == 4)
    option = 'value';
  else
    option = varargin{1};
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

  switch (lower (option))
    case {'value'}
      [eu, F] = sp_eval_msh (u, sp, msh);
      F  = reshape (F, [ndim, npts]);
      eu = squeeze (reshape (eu, [sp.ncomp, npts]));

    case {'divergence'}
      [eu, F] = sp_eval_div_msh (u, sp, msh);
      F  = reshape (F, [ndim, npts]);
      eu = reshape (eu, npts);

    case {'curl'}
      [eu, F] = sp_eval_curl_msh (u, sp, msh);
      F  = reshape (F, [ndim, npts]);
      if (ndim == 2)
        eu = reshape (eu, npts);
      elseif (ndim == 3)
        eu = reshape (eu, [ndim, npts]);
      end

    case {'gradient'}
      [eu, F] = sp_eval_grad_msh (u, sp, msh);
      F  = reshape (F, [ndim, npts]);
      eu = squeeze (reshape (eu, [sp.ncomp, ndim, npts]));

    otherwise
      error ('sp_eval: unknown option to evaluate')
  end

end
