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
% Copyright (C) 2011, 2012, 2014 Rafael Vazquez
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

% Temporary solution, to be fixed using "isprop" after defining the
%  classes with classdef
  knt = cell (ndim, 1);
  if (isfield (struct(space), 'knots'))
    for idim=1:ndim
      knt{idim} = space.knots{idim}(space.degree(idim)+1:end-space.degree(idim));
    end
  elseif (isfield (struct(space), 'sp1'))
    for idim=1:ndim
      knt{idim} = space.sp1.knots{idim}(space.sp1.degree(idim)+1:end-space.sp1.degree(idim));
    end
  else
    for idim=1:ndim; knt{idim} = [0 1]; end
  end
  
  if (iscell (npts))
    pts = npts;
    npts = cellfun (@numel, pts);
  elseif (isvector (npts))
    for idim = 1:ndim
      pts{idim} = linspace (knt{idim}(1), knt{idim}(end), npts(idim));
    end
%     if (ndim == 2)
%       pts = {(linspace (knt{1}(1), knt{1}(end), npts(1))), (linspace (knt{2}(1), knt{2}(end), npts(2)))};
%     elseif (ndim == 3)
%       pts = {(linspace (knt{1}(1), knt{1}(end), npts(1))), (linspace (knt{2}(1), knt{2}(end), npts(2))), (linspace (knt{3}(1), knt{3}(end), npts(3)))};
%     end
  end

  for jj = 1:ndim
    pts{jj} = pts{jj}(:)';
    if (numel (pts{jj}) > 1)
      brk{jj} = [knt{jj}(1), pts{jj}(1:end-1) + diff(pts{jj})/2, knt{jj}(end)];
    else
      brk{jj} = [knt{jj}(1) knt{jj}(end)];
    end
  end

  msh = msh_geopdes (brk, pts, [], geometry, 'boundary', false);
  sp  = space.constructor (msh);

  switch (lower (option))
    case {'value'}
      [eu, F] = sp_eval_msh (u, sp, msh);
      F  = reshape (F, [msh.rdim, npts]);
      eu = squeeze (reshape (eu, [sp.ncomp, npts]));

    case {'divergence'}
      [eu, F] = sp_eval_div_msh (u, sp, msh);
      F  = reshape (F, [msh.rdim, npts]);
      eu = reshape (eu, npts);

    case {'curl'}
      [eu, F] = sp_eval_curl_msh (u, sp, msh);
      F  = reshape (F, [msh.rdim, npts]);
      if (ndim == 2 && msh.rdim == 2)
        eu = reshape (eu, npts);
      elseif (ndim == 3 && msh.rdim == 3)
        eu = reshape (eu, [ndim, npts]);
      else
        error ('sp_eval: the evaluation of the curl for 3d surfaces is not implemented yet')
      end

    case {'gradient'}
      [eu, F] = sp_eval_grad_msh (u, sp, msh);
      F  = reshape (F, [msh.rdim, npts]);
      eu = squeeze (reshape (eu, [sp.ncomp, msh.rdim, npts]));

    otherwise
      error ('sp_eval: unknown option to evaluate')
  end

end
