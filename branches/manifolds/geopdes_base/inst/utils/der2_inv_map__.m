% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
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

function [uxx, uxy, uyy, vxx, vxy, vyy, ux, uy, vx, vy] = der2_inv_map__ (xu, xv, yu, yv, xuu, xuv, xvv, yuu, yuv, yvv)

det = -(xv.*yu - xu.*yv);
det2 = det.^2;

ux =  yv./det;
uy = -xv./det;
vx = -yu./det;
vy =  xu./det;

uxx = (yu.*((uy.*yuv + ux.*xuv).*yv - yu.*(uy.*yvv + ux.*xvv)) - ...
       yv.*((uy.*yuu + ux.*xuu).*yv - yu.*(uy.*yuv + ux.*xuv))) ./ det2;
uxy = (yu.*(xu.*(uy.*yvv + ux.*xvv) - xv.*(uy.*yuv + ux.*xuv)) - ...
       (xu.*(uy.*yuv + ux.*xuv) - xv.*(uy.*yuu + ux.*xuu)).*yv) ./ det2;
uyy = (xv.*(xu.*(uy.*yuv + ux.*xuv) - xv.*(uy.*yuu + ux.*xuu)) - ...
       xu.*(xu.*(uy.*yvv + ux.*xvv) - xv.*(uy.*yuv + ux.*xuv))) ./ det2;

vxx = (yu.*((vy.*yuv + vx.*xuv).*yv - yu.*(vy.*yvv + vx.*xvv)) - ...
       yv.*((vy.*yuu + vx.*xuu).*yv - yu.*(vy.*yuv + vx.*xuv))) ./ det2;
vxy = (yu.*(xu.*(vy.*yvv + vx.*xvv) - xv.*(vy.*yuv + vx.*xuv)) - ...
       (xu.*(vy.*yuv + vx.*xuv) - xv.*(vy.*yuu + vx.*xuu)).*yv) ./ det2;
vyy= (xv.*(xu.*(vy.*yuv + vx.*xuv) - xv.*(vy.*yuu + vx.*xuu)) - ...
      xu.*(xu.*(vy.*yvv + vx.*xvv) - xv.*(vy.*yuv + vx.*xuv))) ./ det2;

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
%! u = rand (100, 1);
%! v = rand (100, 1);
%! [uxx, uxy, uyy, vxx, vxy, vyy, ux, uy, vx, vy] = ...
%!     der2_inv_map__ (xufun (u,v), xvfun (u,v), yufun (u,v), 
%! 		    yvfun (u,v), xuufun (u,v), xuvfun (u,v), 
%! 		    xvvfun (u,v), yuufun (u,v), yuvfun (u,v), 
%! 		    yvvfun (u,v));
%! assert (uxfun (u, v), ux, 1e-10)
%! assert (vyfun (u, v), vy, 1e-10)
%! assert (uxxfun (u, v), uxx, 1e-10)
%! assert (vyyfun (u, v), vyy, 1e-10)