function uknt = kntunclamp (knt, deg, k, dim)
% KNTUNCLAMP: Compute the unclamped knot vector starting from an open one.
%
% Calling Sequence:
% 
%   uknt = kntunclamp (knt, deg, k)
%   uknt = kntunclamp (knt, deg, k, dim)
% 
% INPUT:
% 
%   knt	: open knot vector: see kntrefine
%   deg : polynomial degree of the spline space
%   k   : continuity for the unclamping (from 0 up to p-1)
%   dim : dimensions in which to unclamp (all by default).
%
% OUTPUT:
% 
%   uknt: unclamped knot vector, see nrbmak
% 
% Description:
% 
%     Unclamps directly the open knot vector. See nrbunclamp
%    for further information.
% 
%    Copyright (C) 2013, 2014 Rafael Vazquez
%    Copyright (C) 2020, Bernard Kapidani
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

  knt_is_cell = true;
  if (~iscell (knt))
    knt = {knt};
    knt_is_cell = false;
  end
  uknt = knt;
  
  ndim = numel (knt);
  if (nargin < 4)
    dim = 1:ndim;
  end
  

% if (iscell (knt))
  if (numel(k) < ndim)
    k = [k(:).', k(end) * ones(1, ndim-numel(k))];
  end
  
  assert (numel(deg) == ndim, 'degrees and knots must have the same size');
  for idim = dim
      
    U  = knt{idim};
    
    p  = deg(idim);
    n  = numel(U) - p - 1;
    m  = n + p + 1;
    kk = k(idim);

    if (kk >= p)
      warning ('Taking the maximum k allowed, degree - 1')
      kk = p - 1;
    end

  % Unclamp at left end
    for ii=0:kk
      U(kk-ii+1) = U(kk-ii+2) - (U(n+1-ii) - U(n-ii));
    end

  % Unclamp at right end
    for ii=0:kk
      U(m-kk+ii) = U(m-kk+ii-1) + U(p+ii+1+1) - U(p+ii+1);
    end
    
    uknt{idim} = U;

  end
  
  if (~knt_is_cell)
    uknt = uknt{1};
  end
  
end


