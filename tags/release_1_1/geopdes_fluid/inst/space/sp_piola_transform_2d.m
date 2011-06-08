% SP_PIOLA_TRANSFORM_2D: Map a function space to the physical domain with Piola transform.
%
%     sp = sp_piola_transform_2d (sp, msh)
%
% INPUTS:
%
%     sp:       function space mapped to the physical domain via componentwise mapping.
%     msh:      structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%
% OUTPUT:
%
%    sp: structure representing the discrete function space, with the following fields:
%
%        FIELD_NAME      (SIZE)                              DESCRIPTION
%        ndof            (scalar)                               total number of degrees of freedom
%        nsh_max         (scalar)                               maximum number of shape functions per element
%        nsh             (1 x msh.nel vector)                   actual number of shape functions per each element
%        connectivity    (nsh_max x msh.nel vector)             indices of basis functions that do not vanish in each element
%        shape_functions (2 x msh.nqn x nsh_max x msh.nel)      vector-valued basis functions evaluated at each quadrature node in each element
%        shape_function_gradients
%                        (2 x 2 x msh.nqn x nsh_max x msh.nel)  basis functions gradients evaluated at each quadrature node in each element
%        boundary        (1 x 4 struct array)              struct array representing the space of traces of basis functions on each edge
%        spfun           (function handle)                 function to evaluate an element of the discrete function space, given the Fourier coefficients and a set of points in the parametric space
%
%   For more details, see the documentation
%
% Copyright (C) 2010 Carlo de Falco
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.
% Author: Carlo de Falco <cdf AT users.sourceforge.net>
% Created: 2010-07-21

function spv = sp_piola_transform_2d (spv, msh)
spv = do_sp_piola_transform_2d__ (spv, msh);
if (isfield (msh, 'boundary'))
  for iside = 1:4
    spv.boundary(iside) = do_sp_piola_transform_2d__ (spv.boundary(iside), msh.boundary(iside));
  end
end
end

function spv = do_sp_piola_transform_2d__ (spv, msh)

shape_functions = geopdes_prod__ (msh.geo_map_jac, spv.shape_functions);

xu = squeeze (msh.geo_map_jac(1, 1, :, :));
xv = squeeze (msh.geo_map_jac(1, 2, :, :));
yu = squeeze (msh.geo_map_jac(2, 1, :, :));
yv = squeeze (msh.geo_map_jac(2, 2, :, :));

jacdet = (xu .* yv - xv .* yu);

if (isfield (spv, 'shape_function_gradients'))
  
  if (isfield (msh, 'geo_map_der2'))
    xuu = reshape (msh.geo_map_der2(1, 1, 1, :, :), [], msh.nel);
    yuu = reshape (msh.geo_map_der2(2, 1, 1, :, :), [], msh.nel);
    xuv = reshape (msh.geo_map_der2(1, 1, 2, :, :), [], msh.nel);
    yuv = reshape (msh.geo_map_der2(2, 2, 1, :, :), [], msh.nel);
    xvv = reshape (msh.geo_map_der2(1, 2, 2, :, :), [], msh.nel);
    yvv = reshape (msh.geo_map_der2(2, 2, 2, :, :), [], msh.nel);
    
    wh  = squeeze (spv.shape_functions(1, :, :, :));
    zh  = squeeze (spv.shape_functions(2, :, :, :));
    whu = squeeze (spv.shape_function_gradients(1, 1, :, :, :));
    whv = squeeze (spv.shape_function_gradients(1, 2, :, :, :));
    zhu = squeeze (spv.shape_function_gradients(2, 1, :, :, :));
    zhv = squeeze (spv.shape_function_gradients(2, 2, :, :, :));
    
    for ii=1:spv.nsh_max      
      [wx, wy, zx, zy] = piola_transform_grad_2d__ (xuu, yuu, xuv, yuv, xvv, ... 
				                    yvv, xu, xv, yu, yv, ...
				                    squeeze (wh(:,ii,:)), squeeze (zh(:,ii,:)), ...
				                    squeeze (whu(:,ii,:)), squeeze (whv(:,ii,:)), ...
				                    squeeze (zhu(:,ii,:)), squeeze (zhv(:,ii,:)), jacdet);    
      
      spv.shape_function_gradients(1, 1, :, ii, :) = wx;
      spv.shape_function_gradients(1, 2, :, ii, :) = wy;
      spv.shape_function_gradients(2, 1, :, ii, :) = zx;
      spv.shape_function_gradients(2, 2, :, ii, :) = zy;
    end %endfor

    spv.shape_function_divs = squeeze (spv.shape_function_gradients(1,1,:,:,:) + spv.shape_function_gradients(2,2,:,:,:)); 

  elseif (isfield (spv, 'shape_function_divs')) 

    for ii=1:spv.nsh_max
      spv.shape_function_divs(:,ii,:) = reshape (spv.shape_function_divs(:,ii,:), size (jacdet))./jacdet;
    end
    spv = rmfield (spv, 'shape_function_gradients');
  
  else
    
    spv = rmfield (spv, 'shape_function_gradients');
  
  end %endif
  
end %endif

for ii=1:spv.nsh_max
  spv.shape_functions(1,:,ii,:) = reshape (shape_functions(1,:,ii,:), size (jacdet))./jacdet;
  spv.shape_functions(2,:,ii,:) = reshape (shape_functions(2,:,ii,:), size (jacdet))./jacdet;
end


end %endfunction

function  [wx, wy, zx, zy] = piola_transform_grad_2d__ (xuu, yuu, xuv, yuv, xvv, yvv, xu, xv, yu, yv, wh, zh, whu, whv, zhu, zhv, det);

detdx = (yu.*(yu.*xvv + xv.*yuv - xu.*yvv) + ...
         yv.*(yv.*xuu + xu.*yuv - xv.*yuu) - 2*yu.*yv.*xuv)./det;
detdy = (xu.*(yv.*xuv + xu.*yvv - yu.*xvv) + ...
         xv.*(xuv.*yu + yuu.*xv - yv.*xuu) - 2*xu.*xv.*yuv)./det;
det2 = det.^2;

wx = (-(xu.*wh + xv.*zh).*detdx + ...
      yv.*(xuu.*wh + xuv.*zh + xu.*whu + xv.*zhu) ...
      -yu.*(xuv.*wh + xvv.*zh + xu.*whv + xv.*zhv))./det2;

wy = (-(xu.*wh + xv.*zh).*detdy + ...
      xu.*(xuv.*wh + xvv.*zh + xu.*whv + xv.*zhv) ...
      -xv.*(xuu.*wh + xuv.*zh + xu.*whu + xv.*zhu))./det2;

zx = (-(yu.*wh + yv.*zh).*detdx + ...
      yv.*(yuu.*wh + yuv.*zh + yu.*whu + yv.*zhu) ...
      -yu.*(yuv.*wh + yvv.*zh + yu.*whv + yv.*zhv))./det2;

zy = (-(yu.*wh + yv.*zh).*detdy + ...
      xu.*(yuv.*wh + yvv.*zh + yu.*whv + yv.*zhv) ...
      -xv.*(yuu.*wh + yuv.*zh + yu.*whu + yv.*zhu))./det2;

end
