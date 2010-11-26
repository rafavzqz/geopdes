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

function PI = b2nst__ (U, V, dU, dV, msh, space)

  %% FIXME only order 3 supported at the moment
  if ([dU dV] ~= [3 3])
   error ('b2nst__: only order 3 supported at the moment')
  end
  
  mcp = numel (U) - dU - 2;
  ncp = numel (V) - dV - 2;
  UT = U;
  VT = V;
  UT([5 end-4]) = [];
  VT([5 end-4]) = [];

  %% get numbers of shape functions to replace
  ns = 0;
  for jj = [1 2 4 5]
    ns=ns+1; rmshp(ns) = sub2ind ([mcp+1, ncp+1], 3, jj);
    ns=ns+1; rmshp(ns) = sub2ind ([mcp+1, ncp+1], mcp-1, jj);
  end
  
  for jj = ncp+2-[1 2 4 5];
    ns=ns+1; rmshp(ns) = sub2ind ([mcp+1, ncp+1], 3,   jj);
    ns=ns+1; rmshp(ns) = sub2ind ([mcp+1, ncp+1], mcp-1, jj);
  end

  %% get numbers of shape functions to remove
  ns=ns+1; rmshp(ns) = sub2ind ([mcp+1, ncp+1], 3,     3);
  ns=ns+1; rmshp(ns) = sub2ind ([mcp+1, ncp+1], 3,     ncp-1);
  ns=ns+1; rmshp(ns) = sub2ind ([mcp+1, ncp+1], mcp-1, ncp-1);
  ns=ns+1; rmshp(ns) = sub2ind ([mcp+1, ncp+1], mcp-1, 3);

  nrmshp = numel (rmshp);

  %% construct projection matrix
  P = speye (space.ndof);
  %% remove basis functions to be removed
  P(:, rmshp) = [];

  uv = [reshape(msh.quad_nodes(1, :, :),[],1), reshape(msh.quad_nodes(2, :, :),[],1)]';
  ns = 0;
  for jj = [1 2 3 4]
    ns=ns+1; fun(:,ns) = tbasisfun (uv, [3, 3], {U(3+[0:4]), VT(jj+[0:4])});
    ns=ns+1; fun(:,ns) = tbasisfun (uv, [3, 3], {U(3+[0:4]), VT(end+1-jj-[4:-1:0])});
    ns=ns+1; fun(:,ns) = tbasisfun (uv, [3, 3], {U(end-2-[4:-1:0]), VT(jj+[0:4])});
    ns=ns+1; fun(:,ns) = tbasisfun (uv, [3, 3], {U(end-2-[4:-1:0]), VT(end+1-jj-[4:-1:0])});
  end
  
  M = op_u_v (space, space, msh, ones (size (msh.quad_weights)));

  for jj = 1:ns
    PI(:, jj) = M \ op_f_v (space, msh, reshape (fun(:,jj), size (msh.quad_weights)));
  end

  PI = [P PI];
  
end
