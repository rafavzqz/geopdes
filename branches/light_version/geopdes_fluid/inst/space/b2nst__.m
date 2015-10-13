% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
%
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2011 Rafael Vazquez
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

function PI = b2nst__ (space, knots, degree, msh)

  U = knots{1};
  V = knots{2};
  dU = degree(1);
  dV = degree(2);

  %% FIXME only odd degree supported at the moment
  if (mod (dU, 2) == 0 || mod (dV, 2) == 0)
   error ('b2nst__: only odd degree supported at the moment')
  end
  
  mcp = numel (U) - dU - 1;
  ncp = numel (V) - dV - 1;
  UT = U;
  VT = V;
  UT([dU+2 end-dV-1]) = [];
  VT([dV+2 end-dV-1]) = [];

  %% get numbers of shape functions to replace
  ucenter = (dU+1)/2+1;
  vcenter = (dV+1)/2+1;
  aux = setdiff (1:dV+2, vcenter);
  ns = 0;
  for jj = aux
    ns=ns+1; rmshp(ns) = sub2ind ([mcp, ncp], ucenter, jj);
    ns=ns+1; rmshp(ns) = sub2ind ([mcp, ncp], mcp+1-ucenter, jj);
  end
  
  for jj = ncp+1-aux;
    ns=ns+1; rmshp(ns) = sub2ind ([mcp, ncp], ucenter, jj);
    ns=ns+1; rmshp(ns) = sub2ind ([mcp, ncp], mcp+1-ucenter, jj);
  end

  %% get numbers of shape functions to remove
  ns=ns+1; rmshp(ns) = sub2ind ([mcp, ncp], ucenter, vcenter);
  ns=ns+1; rmshp(ns) = sub2ind ([mcp, ncp], ucenter, ncp+1-vcenter);
  ns=ns+1; rmshp(ns) = sub2ind ([mcp, ncp], mcp+1-ucenter, ncp+1-vcenter);
  ns=ns+1; rmshp(ns) = sub2ind ([mcp, ncp], mcp+1-ucenter, vcenter);

  nrmshp = numel (rmshp);

  %% construct projection matrix
  P = speye (space.ndof);
  %% remove basis functions to be removed
  P(:, rmshp) = [];

 % uv = [reshape(msh.quad_nodes(1, :, :),[],1), reshape(msh.quad_nodes(2, :, :),[],1)]';
  ns = 0;
  for jj = 1:dV+1
    ns=ns+1; fun{ns} = @(pts) tbasisfun (pts, [dU, dV], ...
                 {U(ucenter+[0:dU+1]), VT(jj+[0:dV+1])});
    ns=ns+1; fun{ns} = @(pts) tbasisfun (pts, [dU, dV], ...
                 {U(ucenter+[0:dU+1]), VT(end+1-jj-[dV+1:-1:0])});
    ns=ns+1; fun{ns} = @(pts) tbasisfun (pts, [dU, dV], ...
                 {U(end-ucenter+1-[dU+1:-1:0]), VT(jj+[0:dV+1])});
    ns=ns+1; fun{ns} = @(pts) tbasisfun (pts, [dU, dV], ...
                 {U(end-ucenter+1-[dU+1:-1:0]), VT(end+1-jj-[dV+1:-1:0])});
  end

% Computation of the projection
%  We cannot use the op_f_v_tp operator, because our function handles fun{jj}
%  are for points in the parametric domain, not in the physical domain
  M = op_u_v_tp (space, space, msh, @(x, y) ones (size(x)));

  rhs = zeros (space.ndof, ns);
  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col, 'gradient', false);
    uv = [reshape(msh_col.quad_nodes(1,:,:),[],1), reshape(msh_col.quad_nodes(2,:,:),[],1)]';
    for jj = 1:ns
      rhs(:,jj) = rhs(:,jj) + op_f_v (sp_col, msh_col, ...
                               reshape (fun{jj}(uv), msh_col.nqn, msh_col.nel));
    end
  end

  for jj = 1:ns
    PI(:, jj) = M \ rhs(:,jj);
%    PI(:, jj) = M \ op_f_v_tp (space, msh, reshape (fun(:,jj), size (msh.quad_weights)));
  end

  PI = [P PI];
  
end
