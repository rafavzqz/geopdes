% SP_EVAL_3D:  Evaluate a 3d function at a given set of points.
%
%   [eu, F] = sp_eval_3d (u, space, geometry, npts);
%   [eu, F] = sp_eval_3d (u, space, msh);
%   [eu, F] = sp_eval_3d (u, space, msh, opt);
%
% INPUT:
%
%     u:         vector of dof weights
%     space:     structure representing the space of discrete functions (see sp_bspline_3d_phys)
%     geometry:  geometry structure (see geo_load)
%     pts:       coordinates of points along each parametric direction
%     msh:       msh structure
%     opt:       if the option 'recompute' is added, the values of the shape functions in the space structure are recomputed. By default the option is off.
%
% OUTPUT:
%
%     eu: the function evaluated in the given points
%     F:  grid points in the physical domain, that is, the mapped points
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
%
% Author: Carlo de Falco <cdf AT users.sourceforge.net>

function [eu, F] = sp_eval_3d (u, space, varargin);

  if ((length (varargin) == 2) && ~ischar(varargin{2}))
    geometry = varargin{1};
    pts      = varargin{2};

    for jj = 1:3
      pts{jj} = pts{jj}(:)';
      if (numel (pts{jj}) > 1)
        brk{jj} = [0, pts{jj}(1:end-1) + diff(pts{jj})/2, 1];
      else
        brk{jj} = [0 1];
      end
    end
    
    warn = warning ('query');
    warning off
    msh = msh_3d_tensor_product ({brk{1}, brk{2}, brk{3}}, pts, [], 'no boundary');
    warning(warn);
    msh = msh_push_forward_3d (msh, geometry);
    sp  = feval (space.spfun, msh);
  elseif (length (varargin) == 1)
    msh   = varargin{1};
    sp  = space;
  elseif (ischar(varargin{2}) && strcmpi(varargin{2}, 'recompute'))
    msh   = varargin{1};
    sp  = feval (space.spfun, msh);
  else
    error ('sp_eval_3d: wrong input parameters')
  end

  sp  = feval (space.spfun, msh);

  uc = zeros (size (sp.connectivity));
  c  = sp.connectivity~=0;  uc(c) = u(sp.connectivity(c)); 
  weight = repmat (reshape (uc, [1, 1, sp.nsh_max, msh.nel]), [space.ncomp, msh.nqn, 1, 1]);
  if (space.ncomp == 1)
    saux = size(sp.shape_functions); 
    sp.shape_functions = reshape (sp.shape_functions, [1 saux]);
  end

  eu = squeeze (sum (weight .* sp.shape_functions, 3));
  F  = msh.geo_map;
  if ((length (varargin) == 2) && ~ischar(varargin{2}))
    eu = squeeze (reshape (eu, space.ncomp, numel (pts{1}), numel (pts{2}), numel (pts{3})));
    F  = reshape (F, 3, numel (pts{1}), numel (pts{2}), numel (pts{3}));
  else
    eu = squeeze (reshape (eu, [space.ncomp, size(msh.quad_weights)]));
  end

end


