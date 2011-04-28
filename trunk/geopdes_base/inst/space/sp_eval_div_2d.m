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
    msh = msh_2d_tensor_product ({brk{1}, brk{2}}, pts, [], 'no boundary');
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
      divs = reshape (sp.shape_function_gradients (1, 1, :, :, :) + sp.shape_function_gradients (2, 2, :, :, :), [msh.nqn, sp.nsh_max, msh.nel]);
    end
  else
    divs = sp.shape_function_divs;
  end

  uc = zeros (size (sp.connectivity));
  uc(sp.connectivity~=0) = u(sp.connectivity(sp.connectivity~=0));
  weight = (repmat (reshape (uc, [1, sp.nsh_max, msh.nel]), [msh.nqn, 1, 1]));
 

  div = squeeze (sum (weight .* divs, 2));
  F   = msh.geo_map;
  if ((length (varargin) == 2) && ~ischar(varargin{2}))
    div = squeeze (reshape (div, numel (pts{1}), numel (pts{2})));
    F  = reshape (F, 2, numel (pts{1}), numel (pts{2}));
  else
    div = squeeze (reshape (div, size (msh.quad_weights)));
  end

end

%!test
%! degree       = [4 4];     
%! regularity   = [3 3];   
%! nsub         = [10 10];  
%! nquad        = [4 4];
%! 
%! 
%! geometry = geo_load ('geo_ring.txt');
%! nurbs       = geometry.nurbs;
%! degelev    = max (degree - (nurbs.order-1), 0);
%! nurbs       = nrbdegelev (nurbs, degelev);
%! [rknots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);
%! 
%! nurbs       = nrbkntins (nurbs, nknots);
%! geometry    = geo_load (nurbs);
%! 
%! rule        = msh_gauss_nodes (nquad);
%! [qn, qw]    = msh_set_quad_nodes (geometry.nurbs.knots, rule);
%! msh         = msh_2d_tensor_product (geometry.nurbs.knots, qn, qw);
%! msh         = msh_push_forward_2d (msh, geometry);
%! 
%! sp_scalar = sp_nurbs_2d_phys (nurbs, msh);
%! sp = sp_scalar_to_vector_2d (sp_scalar, sp_scalar, msh, 'divergence', true);
%! 
%! fx  = @(x, y)  reshape(x, [1, size(x)]);
%! fy  = @(x, y)  reshape(-y, [1, size(x)]);
%! f   = @(x, y)  cat (1, fx(x, y), fy(x, y));
%! 
%! [x, y] = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)));
%! mat     = op_u_v (sp, sp, msh, ones (size (x))); 
%! rhs     = op_f_v (sp, msh, f (x, y));
%! 
%! u = mat \ rhs;
%! 
%! vtk_pts = {linspace(0, 1, 31), linspace(0, 1, 31)};
%! [div, F] = sp_eval_div_2d (u, sp, geometry, vtk_pts);
%!
%! assert (norm (div(:), inf), 0, 1e-10)


%!demo
%! degree       = [4 4];     
%! regularity   = [3 3];   
%! nsub         = [10 10];  
%! nquad        = [4 4];
%! 
%! 
%! geometry = geo_load ('geo_ring.txt');
%! nurbs       = geometry.nurbs;
%! degelev    = max (degree - (nurbs.order-1), 0);
%! nurbs       = nrbdegelev (nurbs, degelev);
%! [rknots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);
%! 
%! nurbs       = nrbkntins (nurbs, nknots);
%! geometry    = geo_load (nurbs);
%! 
%! rule        = msh_gauss_nodes (nquad);
%! [qn, qw]    = msh_set_quad_nodes (geometry.nurbs.knots, rule);
%! msh         = msh_2d_tensor_product (geometry.nurbs.knots, qn, qw);
%! msh         = msh_push_forward_2d (msh, geometry);
%! 
%! sp_scalar = sp_nurbs_2d_phys (nurbs, msh);
%! sp = sp_scalar_to_vector_2d (sp_scalar, sp_scalar, msh, 'divergence', true);
%! 
%! fx  = @(x, y)  reshape(x, [1, size(x)]);
%! fy  = @(x, y)  reshape(-y, [1, size(x)]);
%! f   = @(x, y)  cat (1, fx(x, y), fy(x, y));
%! 
%! [x, y] = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)));
%! mat     = op_u_v (sp, sp, msh, ones (size (x))); 
%! rhs     = op_f_v (sp, msh, f (x, y));
%! 
%! u = mat \ rhs;
%! 
%! vtk_pts = {linspace(0, 1, 31), linspace(0, 1, 31)};
%! [div, F] = sp_eval_div_2d (u, sp, geometry, vtk_pts);
%! [X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
%! 
%! 
%! surf (X, Y, div);
%! title ('div \{x, -y\}^T')

