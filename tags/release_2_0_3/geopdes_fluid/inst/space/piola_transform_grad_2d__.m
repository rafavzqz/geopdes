% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
%
% Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
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

function  [wx, wy, zx, zy] = piola_transform_grad_2d__ (xuu, yuu, xuv, yuv, xvv, yvv, xu, xv, yu, yv, wh, zh, whu, whv, zhu, zhv, det);

detdx = (yu.*(yu.*xvv + xv.*yuv - xu.*yvv) + ...
         yv.*(yv.*xuu + xu.*yuv - xv.*yuu) - 2*yu.*yv.*xuv)./det;
detdy = (xu.*(yv.*xuv + xu.*yvv - yu.*xvv) + ...
         xv.*(xuv.*yu + yuu.*xv - yv.*xuu) - 2*xu.*xv.*yuv)./det;
det2 = det.^2;

wx = (-(xu.*wh + xv.*zh).*detdx +...
            yv.*(xuu.*wh + xuv.*zh + xu.*whu + xv.*zhu)...
           -yu.*(xuv.*wh + xvv.*zh + xu.*whv + xv.*zhv))./det2;

wy = (-(xu.*wh + xv.*zh).*detdy +...
            xu.*(xuv.*wh + xvv.*zh + xu.*whv + xv.*zhv)...
           -xv.*(xuu.*wh + xuv.*zh + xu.*whu + xv.*zhu))./det2;

zx = (-(yu.*wh + yv.*zh).*detdx +...
            yv.*(yuu.*wh + yuv.*zh + yu.*whu + yv.*zhu)...
           -yu.*(yuv.*wh + yvv.*zh + yu.*whv + yv.*zhv))./det2;

zy = (-(yu.*wh + yv.*zh).*detdy +...
            xu.*(yuv.*wh + yvv.*zh + yu.*whv + yv.*zhv)...
           -xv.*(yuu.*wh + yuv.*zh + yu.*whu + yv.*zhu))./det2;

end
