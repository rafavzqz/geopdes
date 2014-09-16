% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
%
% Copyright (C) 2009 Rafael Vazquez
% Copyright (C) 2014 Elena Bulgarello, Carlo de Falco, Sara Frizziero
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

function [v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z] = piola_transform_grad_3d__ ...
    (xuu, yuu, zuu, xuv, yuv, zuv, xuw, yuw, zuw, xvv, yvv, zvv, ...
     xvw, yvw, zvw, xww, yww, zww, xu, xv, xw, yu, yv, yw, zu, zv, zw, ...
     u1h, u2h, u3h, u1hu, u1hv, u1hw, u2hu, u2hv, u2hw, u3hu, u3hv, u3hw);

det = xu.*yv.*zw + xv.*yw.*zu + xw.*yu.*zv - ...
      (xw.*yv.*zu + xu.*yw.*zv + xv.*yu.*zw);

ux = (yv.*zw - yw.*zv)./det;
uy = (xw.*zv - xv.*zw)./det;
uz = (xv.*yw - xw.*yv)./det;
vx = (yw.*zu - yu.*zw)./det;
vy = (xu.*zw - xw.*zu)./det;
vz = (xw.*yu - xu.*yw)./det;
wx = (yu.*zv - yv.*zu)./det;
wy = (xv.*zu - xu.*zv)./det;
wz = (xu.*yv - xv.*yu)./det;

xdux = xuu.*ux + xuv.*vx + xuw.*wx;
xdvx = xuv.*ux + xvv.*vx + xvw.*wx;
xdwx = xuw.*ux + xvw.*vx + xww.*wx;
ydux = yuu.*ux + yuv.*vx + yuw.*wx;
ydvx = yuv.*ux + yvv.*vx + yvw.*wx;
ydwx = yuw.*ux + yvw.*vx + yww.*wx;
zdux = zuu.*ux + zuv.*vx + zuw.*wx;
zdvx = zuv.*ux + zvv.*vx + zvw.*wx;
zdwx = zuw.*ux + zvw.*vx + zww.*wx;

detdx = xdux.*(yv.*zw - yw.*zv) + xdvx.*(yw.*zu - yu.*zw) ...
        + xdwx.*(yu.*zv - yv.*zu) + ydux.*(xw.*zv - xv.*zw) ...
        + ydvx.*(xu.*zw - xw.*zu) + ydwx.*(xv.*zu - xu.*zv) ...
        + zdux.*(xv.*yw - xw.*yv) + zdvx.*(xw.*yu - xu.*yw) ...
        + zdwx.*(xu.*yv - xv.*yu);

xduy = xuu.*uy + xuv.*vy + xuw.*wy;
xdvy = xuv.*uy + xvv.*vy + xvw.*wy;
xdwy = xuw.*uy + xvw.*vy + xww.*wy;
yduy = yuu.*uy + yuv.*vy + yuw.*wy;
ydvy = yuv.*uy + yvv.*vy + yvw.*wy;
ydwy = yuw.*uy + yvw.*vy + yww.*wy;
zduy = zuu.*uy + zuv.*vy + zuw.*wy;
zdvy = zuv.*uy + zvv.*vy + zvw.*wy;
zdwy = zuw.*uy + zvw.*vy + zww.*wy;

detdy = xduy.*(yv.*zw - yw.*zv) + xdvy.*(yw.*zu - yu.*zw) ...
        + xdwy.*(yu.*zv - yv.*zu) + yduy.*(xw.*zv - xv.*zw) ...
        + ydvy.*(xu.*zw - xw.*zu) + ydwy.*(xv.*zu - xu.*zv) ...
        + zduy.*(xv.*yw - xw.*yv) + zdvy.*(xw.*yu - xu.*yw) ...
        + zdwy.*(xu.*yv - xv.*yu);

xduz = xuu.*uz + xuv.*vz + xuw.*wz;
xdvz = xuv.*uz + xvv.*vz + xvw.*wz;
xdwz = xuw.*uz + xvw.*vz + xww.*wz;
yduz = yuu.*uz + yuv.*vz + yuw.*wz;
ydvz = yuv.*uz + yvv.*vz + yvw.*wz;
ydwz = yuw.*uz + yvw.*vz + yww.*wz;
zduz = zuu.*uz + zuv.*vz + zuw.*wz;
zdvz = zuv.*uz + zvv.*vz + zvw.*wz;
zdwz = zuw.*uz + zvw.*vz + zww.*wz;

detdz = xduz.*(yv.*zw - yw.*zv) + xdvz.*(yw.*zu - yu.*zw) ...
        + xdwz.*(yu.*zv - yv.*zu) + yduz.*(xw.*zv - xv.*zw) ...
        + ydvz.*(xu.*zw - xw.*zu) + ydwz.*(xv.*zu - xu.*zv) ...
        + zduz.*(xv.*yw - xw.*yv) + zdvz.*(xw.*yu - xu.*yw) ...
        + zdwz.*(xu.*yv - xv.*yu);


u1hx = u1hu.*ux + u1hv.*vx + u1hw.*wx;
u2hx = u2hu.*ux + u2hv.*vx + u2hw.*wx;
u3hx = u3hu.*ux + u3hv.*vx + u3hw.*wx;
u1hy = u1hu.*uy + u1hv.*vy + u1hw.*wy;
u2hy = u2hu.*uy + u2hv.*vy + u2hw.*wy;
u3hy = u3hu.*uy + u3hv.*vy + u3hw.*wy;
u1hz = u1hu.*uz + u1hv.*vz + u1hw.*wz;
u2hz = u2hu.*uz + u2hv.*vz + u2hw.*wz;
u3hz = u3hu.*uz + u3hv.*vz + u3hw.*wz;

v1 = 1./det.*(xu.*u1h + xv.*u2h + xw.*u3h);
v2 = 1./det.*(yu.*u1h + yv.*u2h + yw.*u3h);
v3 = 1./det.*(zu.*u1h + zv.*u2h + zw.*u3h);

v1x = (-detdx.*v1 + ...
(xdux.*u1h + xu.*u1hx + xdvx.*u2h + xv.*u2hx + xdwx.*u3h + xw.*u3hx))./det;
v1y = (-detdy.*v1 + ...
(xduy.*u1h + xu.*u1hy + xdvy.*u2h + xv.*u2hy + xdwy.*u3h + xw.*u3hy))./det;
v1z = (-detdz.*v1 + ...
(xduz.*u1h + xu.*u1hz + xdvz.*u2h + xv.*u2hz + xdwz.*u3h + xw.*u3hz))./det;
v2x = (-detdx.*v2 + ...
(ydux.*u1h + yu.*u1hx + ydvx.*u2h + yv.*u2hx + ydwx.*u3h + yw.*u3hx))./det;
v2y = (-detdy.*v2 + ...
(yduy.*u1h + yu.*u1hy + ydvy.*u2h + yv.*u2hy + ydwy.*u3h + yw.*u3hy))./det;
v2z = (-detdz.*v2 + ...
(yduz.*u1h + yu.*u1hz + ydvz.*u2h + yv.*u2hz + ydwz.*u3h + yw.*u3hz))./det;
v3x = (-detdx.*v3 + ...
(zdux.*u1h + zu.*u1hx + zdvx.*u2h + zv.*u2hx + zdwx.*u3h + zw.*u3hx))./det;
v3y = (-detdy.*v3 + ...
(zduy.*u1h + zu.*u1hy + zdvy.*u2h + zv.*u2hy + zdwy.*u3h + zw.*u3hy))./det;
v3z = (-detdz.*v3 + ...
(zduz.*u1h + zu.*u1hz + zdvz.*u2h + zv.*u2hz + zdwz.*u3h + zw.*u3hz))./det;

end
