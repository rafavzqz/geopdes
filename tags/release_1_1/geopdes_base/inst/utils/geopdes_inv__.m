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

function [Jinv, det] = geopdes_inv__ (v)
  det  = (v(1,1,:,:) .* v(2,2,:,:) - v(2,1,:,:) .* v(1,2,:,:));
  Jinv(1,1,:,:) = v(2,2,:,:)./det;
  Jinv(2,2,:,:) = v(1,1,:,:)./det;
  Jinv(1,2,:,:) = -v(1,2,:,:)./det;
  Jinv(2,1,:,:) = -v(2,1,:,:)./det;
  det = squeeze (det);
end
