% geopdes_inv: the inverse matrix or, for codimension > 0 (rdim > ndim)
% the Moore-Penrose pseudo-inverse.
%
% To be used with msh.geo_map_jac, to compute DF^{-1}.
% In the case of codimension > 0, it computes
%   J = (DF^t * DF)^{-1} * DF^t
%   det = sqrt ( det(DF^t * DF)^{-1} )
%
%  [JinvT, det] = geopdes_invT__ (geo_map_jac)
%
% OUTPUT:
%   Jinv: the computed matrix evaluated at every point. Size (ndim x rdim x nqn x nel)
%   det:  the determinant evaluated at every point. Size (nqn x nel)
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

function [Jinv, det] = geopdes_inv__ (v)

  [JinvT, det] = geopdes_invT__ (v);
  Jinv = permute (JinvT, [2, 1, 3:numel(size(v))]);

end

%!test
%! A = [1 2; 3 4];
%! Ainv = geopdes_inv__ (A);
%! assert (Ainv * A, eye(2), 1e-14)
%!
%!test
%! A = [1 2 3; 4 5 6; 7 8 10];
%! Ainv = geopdes_inv__ (A);
%! assert (Ainv * A, eye(3), 1e-14)
