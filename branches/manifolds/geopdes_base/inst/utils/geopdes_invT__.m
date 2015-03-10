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

function [JinvT, det] = geopdes_invT__ (v)

  vsize = size (v);
  if (numel (vsize) < 3)
    vsize(3) = 1;
  end

  if (vsize(1) == 2 && vsize(2) == 2)
    det = v(1,1,:,:) .* v(2,2,:,:) - v(2,1,:,:) .* v(1,2,:,:);
    
    JinvT(1,1,:,:) = v(2,2,:,:)./det;
    JinvT(2,2,:,:) = v(1,1,:,:)./det;
    JinvT(1,2,:,:) = -v(2,1,:,:)./det;
    JinvT(2,1,:,:) = -v(1,2,:,:)./det;
  elseif (vsize(1) == 3 && vsize(2) == 3)
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

  elseif (vsize(1) == 3 && vsize(2) == 2)
    % G = v^t * v, first fundamental form for v = DF
    for ii = 1:2
      for jj = 1:2
        G(ii,jj,:,:) = sum (v(:,ii,:,:).*v(:,jj,:,:), 1);
      end
    end
    det = sqrt (G(1,1,:,:) .* G(2,2,:,:) - G(2,1,:,:) .* G(1,2,:,:));
    
    Ginv = geopdes_inv__ (G);
    %  v * G^{-1}, which means DF * {DF^t * DF}^(-1)
    JinvT = zeros (size(v));
    for ii = 1:3
      for jj = 1:2
        JinvT(ii,jj,:,:) = sum (reshape (v(ii,:,:,:), vsize(2:end)) .* ...
            reshape (Ginv(:,jj,:,:), [2, vsize(3:end)]), 1);
      end
    end
    
  end
  
  det = squeeze (det);

end
