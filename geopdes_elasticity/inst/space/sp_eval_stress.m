% SP_EVAL_STRESS: Compute the stress of the solution, given by its degrees of freedom, at a given set of points.
%
%   [eu, F] = sp_eval_stress (u, space, geometry, lambda, mu, pts);
%   [eu, F] = sp_eval_stress (u, space, geometry, lambda, mu, npts);
%
% INPUT:
%     
%     u:          vector of dof weights
%     space:      object defining the discrete space (see sp_vector)
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
% Copyright (C) 2011, 2012, 2015 Rafael Vazquez
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
  
  der2 = false;
  if (isa (space, 'sp_vector_piola_transform'))
    der2 = true;
  end

  ndim = numel (npts);

% Temporary solution, to be fixed using "isprop" after defining the
%  classes with classdef
  knt = cell (ndim, 1);
  if (isfield (struct(space), 'knots'))
    for idim=1:ndim
      knt{idim} = space.knots{idim}(space.degree(idim)+1:end-space.degree(idim));
    end
  elseif (isfield (struct(space), 'scalar_spaces'))
    for idim=1:ndim
      knt{idim} = space.scalar_spaces{1}.knots{idim}(space.scalar_spaces{1}.degree(idim)+1:end-space.scalar_spaces{1}.degree(idim));
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
  end

  for jj = 1:ndim
    pts{jj} = pts{jj}(:)';
    if (numel (pts{jj}) > 1)
      brk{jj} = [knt{jj}(1), pts{jj}(1:end-1) + diff(pts{jj})/2, knt{jj}(end)];
    else
      brk{jj} = [knt{jj}(1) knt{jj}(end)];
    end
  end

  msh = msh_cartesian (brk, pts, [], geometry, 'boundary', false, 'der2', der2);
  sp  = space.constructor (msh);

  [stress, F] = sp_eval_stress_msh (u, sp, msh, lambda, mu);
  F  = reshape (F, [ndim, npts]);
  stress = reshape (stress, [sp.ncomp, ndim, npts]);

end
