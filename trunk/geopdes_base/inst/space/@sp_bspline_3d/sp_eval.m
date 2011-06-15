% SP_EVAL: Evaluate a function, given by its degrees of freedom, at a given set of points.
%
%   [eu, F] = sp_eval (u, space, geometry, pts);
%   [eu, F] = sp_eval (u, space, msh);
%   [eu, F] = sp_eval (u, space, msh, opt);
%
% INPUT:
%     
%     u:         vector of dof weights
%     space:     class defining the space (see sp_bspline_3d)
%     geometry:  geometry structure (see geo_load)
%     pts:       coordinates of points along each parametric direction
%     msh:       msh class (see msh_3d)
%     opt:       if the option 'recompute' is added, the values of the shape functions in the space structure are recomputed. By default the option is off.
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

function [eu, F] = sp_eval (u, space, varargin);

  recompute = false;

  if (~rem (length (varargin), 2) == 0)
    method = 'msh';
  else
    method = 'points';
  end

  if (strcmp (method, 'points'))
    geometry = varargin {1};
    pts      = varargin {2};

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
    sp  = sp_bspline_3d (space.knots, space.degree, msh);

  elseif (strcmp (method, 'msh'))
    msh = varargin {1};
    if (recompute)
      sp  = sp_bspline_3d (space.knots, space.degree, msh);
    else
      sp = space;
    end
  else
    error ('sp_eval: wrong number of input parameters')
  end

  F  = msh.geo_map;
  eu = zeros (msh.nqn, msh.nel);

  nel_col = msh.nelv * msh.nelw;
  for iel = 1:msh.nelu
    [sp_col, elem_list] = sp_evaluate_col (sp, msh, iel, 'gradient', false);

    uc_iel = zeros (size (sp_col.connectivity));
    uc_iel(sp_col.connectivity~=0) = ...
          u(sp_col.connectivity(sp_col.connectivity~=0));
    weight = repmat (reshape (uc_iel, [1, sp_col.nsh_max, nel_col]), ...
                                  [msh.nqn, 1, 1]);

    eu(:, elem_list) = squeeze (sum (weight .* sp_col.shape_functions, 2));
  end

  if (strcmp (method, 'points'))
    eu = reshape (eu, numel (pts{1}), numel (pts{2}), numel (pts{3}));
    F  = reshape (F, 3, numel (pts{1}), numel (pts{2}), numel (pts{3}));
  end

end
