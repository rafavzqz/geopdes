% SP_EVAL_DIV_2D: Evaluate the divergence of a 2d filed at a given set of points.
%
%   [div, F] = sp_eval_div_2d (u, space, geometry, npts);
%   [div, F] = sp_eval_div_2d (u, space, msh);
%   [div, F] = sp_eval_div_2d (u, space, msh, opt);
%
% INPUT:
%     
%     u:         vector of dof weights
%     space:     structure representing the space of discrete functions (see sp_bspline_2d_phys)
%     geometry:  geometry structure (see geo_load)
%     pts:       coordinates of points along each parametric direction
%     msh:       msh structure
%     opt:       if the option 'recompute' is added, the values of the shape functions in the space structure are recomputed. By default the option is off.
%
% OUTPUT:
%
%     div: the divergence of the field evaluated in the given points 
%     F:   grid points in the physical domain, that is, the mapped points
% 
% Copyright (C) 2011 Carlo de Falco
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

function [div, F] = sp_eval_div_2d (u, space, varargin);

  if ((length (varargin) == 2) && ~ischar(varargin{2}))
    geometry = varargin {1};
    pts      = varargin {2};
    ku = [0, pts{1}(1:end-1) + diff(pts{1})/2, 1];
    kv = [0, pts{2}(1:end-1) + diff(pts{2})/2, 1];
    
    warn = warning ('query');
    warning off
    msh = msh_2d_tensor_product ({ku, kv}, pts, [], 'no boundary');
    warning(warn);
    msh = msh_push_forward_2d (msh, geometry);
    sp  = feval (space.spfun, msh);
  elseif (length (varargin) == 1)
    msh = varargin {1};
    sp  = space;
  elseif (ischar(varargin{2}) && strcmpi(varargin{2}, 'recompute'))
    msh   = varargin {1};
    sp  = feval (space.spfun, msh);
  else
    error ('sp_eval_div_2d: wrong number of input parameters')
  end

  if (space.ncomp == 1)
    error ('sp_eval_div_2d: field cannot be scalar')
  end

  if (~isfield (sp, 'shape_function_divs'))
    if (~isfield (sp, 'shape_function_gradients'))
      error ('sp_eval_div_2d: the space structure must have either divs or gradients')
    else
      divs = squeeze (sp.shape_function_gradients (1, 1, :, :, :) + sp.shape_function_gradients (2, 2, :, :, :));
    end
  else
    divs = sp.shape_function_divs;
  end

  uc = zeros (size (sp.connectivity));
  uc(sp.connectivity~=0) = u(sp.connectivity(sp.connectivity~=0));
  weight = repmat (reshape (uc, [1, sp.nsh_max, msh.nel]), [msh.nqn, 1, 1]);
 

  div = squeeze (sum (weight .* divs, 2));
  F   = msh.geo_map;
  if ((length (varargin) == 2) && ~ischar(varargin{2}))
    div = squeeze (reshape (div, numel (pts{1}), numel (pts{2})));
    F  = reshape (F, 2, numel (pts{1}), numel (pts{2}));
  else
    div = squeeze (reshape (div, size (msh.quad_weights)));
  end

end
