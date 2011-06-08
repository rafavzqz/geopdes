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

function c = geopdes_prod__ (a, b)
nel = size(b,4);
c = zeros(size(a,1), size(a,3), size(b, 3), nel);
  for ish = 1:size (b, 3)
    for idir = 1:size (a,1)
      for jdir = 1:size(a,2)
        tmp_a = a(idir, jdir, :, :);
        tmp_b = b(jdir, :, ish, :);
	c(idir,:,ish,:) = c(idir,:,ish,:) + ...
	  reshape(tmp_a(:).*tmp_b(:), 1, [], 1, nel);
      end
    end
  end
end
