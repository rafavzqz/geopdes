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

function [JinvT, det] = geopdes_invT__ (v)

  if (size (v,1) == 2)
    det  = (v(1,1,:,:) .* v(2,2,:,:) - v(2,1,:,:) .* v(1,2,:,:));
    JinvT(1,1,:,:) = v(2,2,:,:)./det;
    JinvT(2,2,:,:) = v(1,1,:,:)./det;
    JinvT(1,2,:,:) = -v(2,1,:,:)./det;
    JinvT(2,1,:,:) = -v(1,2,:,:)./det;
    det = squeeze (det);
  elseif (size (v,1) == 3)
    det = v(1,1,:,:) .* (v(2,2,:,:) .* v(3,3,:,:) - v(2,3,:,:) .* v(3,2,:,:))...
        + v(1,2,:,:) .* (v(2,3,:,:) .* v(3,1,:,:) - v(2,1,:,:) .* v(3,3,:,:))...
        + v(1,3,:,:) .* (v(2,1,:,:) .* v(3,2,:,:) - v(2,2,:,:) .* v(3,1,:,:));

    JinvT(1,1,:,:) =  (v(2,2,:,:).*v(3,3,:,:) - v(2,3,:,:).*v(3,2,:,:))./det;
    JinvT(1,2,:,:) = -(v(2,1,:,:).*v(3,3,:,:) - v(2,3,:,:).*v(3,1,:,:))./det;
    JinvT(1,3,:,:) =  (v(2,1,:,:).*v(3,2,:,:) - v(2,2,:,:).*v(3,1,:,:))./det;
    JinvT(2,1,:,:) = -(v(1,2,:,:).*v(3,3,:,:) - v(1,3,:,:).*v(3,2,:,:))./det;
    JinvT(2,2,:,:) =  (v(1,1,:,:).*v(3,3,:,:) - v(1,3,:,:).*v(3,1,:,:))./det;
    JinvT(2,3,:,:) = -(v(1,1,:,:).*v(3,2,:,:) - v(1,2,:,:).*v(3,1,:,:))./det;
    JinvT(3,1,:,:) =  (v(1,2,:,:).*v(2,3,:,:) - v(1,3,:,:).*v(2,2,:,:))./det;
    JinvT(3,2,:,:) = -(v(1,1,:,:).*v(2,3,:,:) - v(1,3,:,:).*v(2,1,:,:))./det;
    JinvT(3,3,:,:) =  (v(1,1,:,:).*v(2,2,:,:) - v(1,2,:,:).*v(2,1,:,:))./det;
    det = squeeze (det);
  end

end
