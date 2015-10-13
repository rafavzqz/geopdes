% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
%
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2015 Rafael Vazquez
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

function [Jinv, geo_inv_der2] = geopdes_inv_der2__ (geo_map_jac, geo_map_der2)

vsize = size (geo_map_jac);
vsize(end+1:4) = 1;

rdim = vsize(1);
ndim = vsize(2);

Jinv = geopdes_inv__ (geo_map_jac);
Jinv = reshape (Jinv, [ndim, rdim, vsize(3:end)]);
Jsize = size (Jinv);
Jsize(end+1:4) = 1;

Jinv1 = reshape (Jinv, [1, Jsize]);
Jinv2 = permute (Jinv1, [1 3 2 4:numel(Jsize)+1]);
Jinv3 = permute (Jinv1, [3 2 1 4:numel(Jsize)+1]);

geo_inv_der2 = zeros ([ndim, rdim, rdim, vsize(3:end)]);

for alp = 1:ndim
  for idim = 1:rdim
    for jdim = 1:rdim
      aux1 = sum (bsxfun (@times, geo_map_der2, Jinv1(:,:,idim,:,:)), 2);
      aux2 = sum (bsxfun (@times, aux1, Jinv2(:,jdim,:,:,:)), 3);
      aux3 = sum (bsxfun (@times, aux2, Jinv3(:,alp,:,:,:)), 1);
      
      geo_inv_der2(alp,idim,jdim,:,:) = -reshape (aux3, vsize(3:end));
    end
  end
end

end

%!test 
%!
%! %% x = u^2 + u + 1 
%! %% y = v^2 + v + 1 
%!
%! xufun = @(u, v) (2*u + 1);
%! xvfun = @(u, v) (0*u);
%! yufun = @(u, v) (0*u);
%! yvfun = @(u, v) (2*v + 1);
%!
%! xuufun = @(u, v) (2+0*u);
%! xuvfun = @(u, v) (0*u);
%! xvvfun = @(u, v) (0*u);
%!
%! yuufun = @(u, v) (0*u);
%! yuvfun = @(u, v) (0*u);
%! yvvfun = @(u, v) (2+0*v);
%!
%! uxfun = @(u, v) (1./(2*u+1));
%! vyfun = @(u, v) (1./(2*v+1));
%! uxxfun = @(u, v) (-2./(2*u+1).^3);
%! vyyfun = @(u, v) (-2./(2*v+1).^3);
%!
%! u = linspace (0,1, 10);
%! v = linspace (0,1, 10);
%! %u = rand (1, 10);
%! %v = rand (1, 10);
%! geo_map_jac(1,1,:) = xufun(u,v);
%! geo_map_jac(1,2,:) = xvfun(u,v);
%! geo_map_jac(2,1,:) = yufun(u,v);
%! geo_map_jac(2,2,:) = yvfun(u,v);
%!
%! geo_map_der2(1,1,1,:) = xuufun(u,v);
%! geo_map_der2(1,1,2,:) = xuvfun(u,v);
%! geo_map_der2(1,2,1,:) = xuvfun(u,v);
%! geo_map_der2(1,2,2,:) = xvvfun(u,v);
%! geo_map_der2(2,1,1,:) = yuufun(u,v);
%! geo_map_der2(2,1,2,:) = yuvfun(u,v);
%! geo_map_der2(2,2,1,:) = yuvfun(u,v);
%! geo_map_der2(2,2,2,:) = yvvfun(u,v);
%!
%! [Jinv, geo_inv_der2] = geopdes_inv_der2__ (geo_map_jac, geo_map_der2);
%! assert (uxfun (u, v), reshape (Jinv(1,1,:), 1, numel(u)), 1e-14)
%! assert (vyfun (u, v), reshape (Jinv(2,2,:), 1, numel(u)), 1e-14)
%! assert (uxxfun (u, v), reshape (geo_inv_der2(1,1,1,:), 1, numel(u)), 1e-14)
%! assert (vyyfun (u, v), reshape (geo_inv_der2(2,2,2,:), 1, numel(u)), 1e-14)