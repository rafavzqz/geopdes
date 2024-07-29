% SP_H2_ERROR: Evaluate the error in H^2 norm, H^1 and L^2 norms.
%
%   [errh2, errh1, errl2, errh2s, errh1s] = sp_h2_error (sp, msh, u, uex, graduex, hessuex)
%
% INPUT:
%
%   space: object defining the space of discrete functions (see sp_multipatch_C1)
%   msh:   object defining the domain partition and the quadrature rule (see msh_multipatch)
%   u:     vector of dof weights
%   uex:   function handle to evaluate the exact solution
%   graduex: function handle to evaluate the gradient of the exact solution
%   hessuex: function handle to evaluate the hessian of the exact solution
%
% OUTPUT:
%
%   errh2: error in H^2 norm
%   errh1: error in H^1 norm
%   errl2: error in L^2 norm
%   errh2s: error in H^2 seminorm
%   errh1s: error in H^1 seminorm
%
% Copyright (C) 2015, 2017, 2022 Rafael Vazquez
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

function [errh2, errh1, errl2, errh2s, errh1s] = sp_h2_error (space, msh, u, uex, graduex, hessuex)

  if (msh.rdim == 3)
    error ('The H2 error is only computed for planar surfaces') 
  end
  if (space.npatch ~= msh.npatch)
    error ('The number of patches does not coincide') 
  end

  for iptc = 1:msh.npatch
%     if (isempty (space.dofs_ornt))
    [Cpatch, Cpatch_cols] = sp_compute_Cpatch (space, iptc);
    u_ptc = Cpatch * u(Cpatch_cols);
%     else
%       u_ptc = u(space.gnum{iptc}) .* space.dofs_ornt{iptc}.';
%     end
    [error_h2(iptc), error_h1(iptc), error_l2(iptc), error_h2s(iptc), error_h1s(iptc)] = ...
        sp_h2_error (space.sp_patch{iptc}, msh.msh_patch{iptc}, u_ptc, uex, graduex, hessuex);
  end
  errl2 = sqrt (sum (error_l2 .* error_l2));
  errh1 = sqrt (sum (error_h1 .* error_h1));
  errh1s = sqrt (sum (error_h1s .* error_h1s));
  errh2 = sqrt (sum (error_h2 .* error_h2));
  errh2s = sqrt (sum (error_h2s .* error_h2s));

end
