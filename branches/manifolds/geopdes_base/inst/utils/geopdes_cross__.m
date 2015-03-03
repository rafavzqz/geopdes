% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
%
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

function w = geopdes_cross__ (u, v)

  w = zeros (size (u));
  w(1,:,:) = u(2,:,:) .* v(3,:,:) - u(3,:,:) .* v(2,:,:);
  w(2,:,:) = u(3,:,:) .* v(1,:,:) - u(1,:,:) .* v(3,:,:);
  w(3,:,:) = u(1,:,:) .* v(2,:,:) - u(2,:,:) .* v(1,:,:);
   
end
