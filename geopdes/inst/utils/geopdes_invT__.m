% geopdes_invT: transpose of the inverse or, for codimension > 0 (rdim > ndim)
% the transpose of the Moore-Penrose pseudo-inverse.
%
% To be used with msh.geo_map_jac, to compute DF^{-T}.
% In the case of codimension > 0, it computes
%   J^T = DF * (DF^t * DF)^{-1}
%   det = sqrt ( det(DF^t * DF)^{-1} )
%
%  [JinvT, det] = geopdes_invT__ (geo_map_jac)
%
% OUTPUT:
%   JinvT: the computed matrix evaluated at every point. Size (rdim x ndim x nqn x nel)
%   det:   the determinatn evaluated at every point. Size (nqn x nel)
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

  elseif (vsize(1) == 1 && vsize(2) == 1)
    det = v;
    JinvT = 1./v;
    
  elseif (vsize(1) > vsize(2))
    % G = v^t * v, first fundamental form for v = DF
    for ii = 1:vsize(2)
      for jj = 1:vsize(2)
        G(ii,jj,:,:) = sum (v(:,ii,:,:).*v(:,jj,:,:), 1);
      end
    end
    [Ginv,det] = geopdes_inv__ (G);
    det = sqrt (det);

    %  v * G^{-1}, which means DF * {DF^t * DF}^(-1)
    JinvT = zeros (size(v));
    for ii = 1:vsize(1)
      for jj = 1:vsize(2)
        JinvT(ii,jj,:,:) = sum (reshape (v(ii,:,:,:), vsize(2:end)) .* ...
            reshape (Ginv(:,jj,:,:), [vsize(2), vsize(3:end)]), 1);
      end
    end
  
  end
  
  det = squeeze (det);

end
