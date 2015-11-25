% SP_EVAL: Compute the value or the derivatives of a function, given by its degrees of freedom, at a given set of points.
%
%   [eu, F] = sp_eval (u, space, geometry, pts, [option], [lambda_lame, mu_lame]);
%   [eu, F] = sp_eval (u, space, geometry, npts, [option], [lambda_lame, mu_lame]);
%
% INPUT:
%     
%     u:           vector of dof weights
%     space:       object defining the discrete space (see sp_vector)
%     geometry:    geometry structure (see geo_load)
%     pts:         cell array with coordinates of points along each parametric direction
%     npts:        number of points along each parametric direction
%     options:     cell array with the fields to plot
%                   accepted options are 'value' (default), 'gradient', 'curl', 'divergence', 'stress'
%     lambda_lame: function handle to the first Lame coefficient (only needed to compute 'stress')
%     mu_lame:     function handle for the second Lame coefficient (only needed to compute 'stress')
%
% OUTPUT:
%
%     eu: cell-array with the fields evaluated at the given points 
%     F:  grid points in the physical domain, that is, the mapped points
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2012, 2014, 2015 Rafael Vazquez
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

function [eu, F] = sp_eval (u, space, geometry, npts, options, lambda_lame, mu_lame)

  if (nargin < 5)
    options = {'value'};
    lambda_lame = [];
    mu_lame = [];
  elseif (nargin < 7)
    lambda_lame = [];
    mu_lame = [];      
  end
  if (~iscell (options))
    options = {options};
  end
  nopts = numel (options);

  ndim = numel (npts);

% Temporary solution, to be fixed using "isprop" after defining the
%  classes with classdef
  knt = cell (ndim, 1);
  if (isfield (struct(space), 'knots'))
    for idim=1:ndim
      knt{idim} = space.knots{idim}(space.degree(idim)+1:end-space.degree(idim));
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

  msh = msh_cartesian (brk, pts, [], geometry, 'boundary', false);
  sp  = space.constructor (msh);

  
  value = false; grad = false; laplacian = false;
  curl = false; divergence = false;
  
  for iopt = 1:nopts
    switch (lower (options{iopt}))
      case 'value'
        eu{iopt} = zeros (space.ncomp, msh.nqn, msh.nel);
        eunum{iopt} = {1:space.ncomp, 1:msh.nqn};
        eusize{iopt} = [space.ncomp, npts];
        value = true;

      case 'gradient'
        eu{iopt} = zeros (space.ncomp, msh.rdim, msh.nqn, msh.nel);
        eunum{iopt} = {1:space.ncomp, 1:msh.rdim, 1:msh.nqn};
        eusize{iopt} = [space.ncomp, msh.rdim, npts];
        grad = true;
        
%       case 'laplacian'
%         eu{iopt} = zeros (space.ncomp, msh.nqn, msh.nel);
%         eunum{iopt} = {1:space.ncomp, 1:msh.nqn};
%         eusize{iopt} = [space.ncomp, npts];
%         laplacian = true;
% 
      case 'curl'
        if (msh.ndim == 2 && msh.rdim == 2)
          eu{iopt} = zeros (msh.nqn, msh.nel);
          eunum{iopt} = {1:msh.nqn};
          eusize{iopt} = npts;
        elseif (msh.ndim == 3 && msh.rdim == 3)
          eu{iopt} = zeros (msh.rdim, msh.nqn, msh.nel);
          eunum{iopt} = {1:msh.rdim, 1:msh.nqn};
          eusize{iopt} = [msh.rdim, npts];
        end
        curl = true;
        
      case 'divergence'
        eu{iopt} = zeros (msh.nqn, msh.nel);
        eunum{iopt} = {1:msh.nqn};
        eusize{iopt} = npts;
        divergence = true;

      case 'stress'
        if (nargin < 6)
          error ('sp_eval_msh: Lame coefficients missing')
        end
        eu{iopt} = zeros (space.ncomp, msh.rdim, msh.nqn, msh.nel);
        eunum{iopt} = {1:space.ncomp, 1:msh.rdim, 1:msh.nqn};
        eusize{iopt} = [sp.ncomp, msh.rdim, npts];
        grad = true; divergence = true;
    end
  end

  F = zeros (msh.rdim, msh.nqn, msh.nel);
  
  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (sp, msh_col, 'value', value, 'gradient', grad, ...
                  'curl', curl, 'divergence', divergence);
% 'laplacian', laplacian, 

    eu_aux = sp_eval_msh (u, sp_col, msh_col, options, lambda_lame, mu_lame);
    
    F(:,:,msh_col.elem_list) = msh_col.geo_map;
    for iopt = 1:nopts
      eu{iopt}(eunum{iopt}{:},msh_col.elem_list) = eu_aux{iopt};
    end
  end
  
  F = reshape (F, [msh.rdim, npts]);
  for iopt = 1:nopts
    eu{iopt} = reshape (eu{iopt}, eusize{iopt});
  end

  if (nopts == 1)
    eu = eu{1};
  end
end
