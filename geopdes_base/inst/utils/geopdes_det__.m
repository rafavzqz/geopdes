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
% Author: Carlo de Falco <cdf AT users.sourceforge.net>

function d = geopdes_det__ (v)  

  vsize = size (v);

  if (vsize(1) == 2 && vsize(2) == 2)
    d = v(1,1,:,:) .* v(2,2,:,:) - v(2,1,:,:) .* v(1,2,:,:);
  elseif (vsize(1) == 3 && vsize(2) == 3)
    d = v(1,1,:,:) .* (v(2,2,:,:) .* v(3,3,:,:) - v(2,3,:,:) .* v(3,2,:,:))...
      + v(1,2,:,:) .* (v(2,3,:,:) .* v(3,1,:,:) - v(2,1,:,:) .* v(3,3,:,:))...
      + v(1,3,:,:) .* (v(2,1,:,:) .* v(3,2,:,:) - v(2,2,:,:) .* v(3,1,:,:));
  elseif (vsize(1) == 3 && vsize(2) == 2)
    % G = v^t * v
    G = zeros (2, 2, vsize(3), vsize(4));
    G(1,1,:,:) = v(1,1,:,:).^2 + v(2,1,:,:).^2 + v(3,1,:,:).^2;
    G(2,2,:,:) = v(1,2,:,:).^2 + v(2,2,:,:).^2 + v(3,2,:,:).^2;
    G(1,2,:,:) = v(1,1,:,:).*v(1,2,:,:) + v(2,1,:,:).*v(2,2,:,:) + v(3,1,:,:).v(3,2,:,:);
    G(2,1,:,:) = G(1,2,:,:);
    
    d = G(1,1,:,:) .* G(2,2,:,:) - G(2,1,:,:) .* G(1,2,:,:);
  end

  d = squeeze (d);
  
end
