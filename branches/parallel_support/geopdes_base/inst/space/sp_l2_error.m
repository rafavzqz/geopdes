% SP_L2_ERROR: Evaluate the error in L^2 norm.
%
%   toterr = sp_l2_error (space, msh, u, uex);
%
% INPUT:
%
%     space:    structure representing the space of discrete functions (see sp_bspline_3d_phys)
%     msh:      structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%     u:        vector of dof weights
%     uex:      function handle to evaluate the exact solution
%
% OUTPUT:
%
%     toterr: error in L^2 norm
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

function toterr = sp_l2_error (sp, msh, u, uex);
  
  if (size (msh.quad_nodes, 1) == 2)
    valnum = sp_eval_2d (u, sp, msh);
  elseif (size (msh.quad_nodes, 1) == 3)
    valnum = sp_eval_3d (u, sp, msh);
  end


  w = msh.quad_weights(:) .* msh.jacdet(:);
  for idir = 1:size (msh.geo_map, 1)
    x{idir} = reshape (msh.geo_map(idir,:,:), msh.nqn*msh.nel, 1);
  end

  valex = feval (uex, x{:});
  valnum = reshape (valnum, sp.ncomp, []);
  valex  = reshape (valex, size (valnum));

  err2 = sum ((valnum - valex).^2, 1);
  toterr = sqrt (sum (err2(:).*w));

end
