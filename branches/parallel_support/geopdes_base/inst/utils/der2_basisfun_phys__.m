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

function [bxx, bxy, byy] = der2_basisfun_phys__ (xu, xv, yu, yv, uxx, uxy, uyy, vxx, vxy, vyy, buu, buv, bvv, bu, bv)

det = (-xv .* yu + xu .* yv);

ux =  yv ./ det;
uy = -xv ./ det;
vx = -yu ./ det;
vy =  xu ./ det;

bxx = buu.*ux.^2 + bvv.*vx.^2 + bu.*uxx + bv.*vxx + 2*buv.*ux.*vx;
byy = buu.*uy.^2 + bvv.*vy.^2 + bu.*uyy + bv.*vyy + 2*buv.*uy.*vy;
bxy = buu.*ux.*uy + bvv.*vx.*vy + bu.*uxy + bv.*vxy + buv.*ux.*vy + buv.*uy.*vx;

end