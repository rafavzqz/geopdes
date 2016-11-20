% Compute the product between matrix A and the 3D tensor X,
%  in the direction d (d = 1, 2 or 3)
%
% Copyright (C) 2016 Mattia Tani
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

function Y = tprod__ (A, X, d)

[nx, ny, nz] = size(X);
m = size(A,1);

if (d == 1)
  Y = A*reshape (X, nx, ny*nz);
  Y = reshape (Y, m, ny, nz);
elseif (d == 2)
  Y = A*reshape (permute (X, [2 1 3]), ny, nx*nz);
  Y = permute (reshape (Y, m, nx, nz), [2 1 3]);
elseif (d == 3)
  Y = A*reshape (permute (X, [3 2 1]), nz, nx*ny);
  Y = permute (reshape (Y, m, ny, nx), [3 2 1]);    
end

end
