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

function [uxx, uxy, uyy, vxx, vxy, vyy] = der2_inv_map__ (xu, xv, yu, yv, xuu, xuv, xvv, yuu, yuv, yvv)

uxx = (-(xv.*yu.*yuv) + xuv.*yu.*yv + xv.*yuu.*yv - xuu.*yv.^2)./(xv.*yu - xu.*yv).^2;
uxy = (-(xv.^2.*yuu) + xu.*xv.*yuv - xu.*xuv.*yv + xuu.*xv.*yv)./(xv.*yu - xu.*yv).^2;
uyy = (-(xv.^2.*yuv) + xuv.*xv.*yv - xu.*xvv.*yv + xu.*xv.*yvv)./(xv.*yu - xu.*yv).^2;

end