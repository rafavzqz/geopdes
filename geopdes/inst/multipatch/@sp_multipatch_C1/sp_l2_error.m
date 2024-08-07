% SP_L2_ERROR: Evaluate the error in L^2 norm.
%
%   errl2 = sp_l2_error (space, msh, u, uex)
%
% INPUT:
%
%   space: object defining the space of discrete functions (see sp_multipatch_C1)
%   msh:   object defining the domain partition and the quadrature rule (see msh_multipatch)
%   u:     vector of dof weights
%   uex:   function handle to evaluate the exact solution
%
% OUTPUT:
%
%     errl2:  error in L^2 norm
%
% Copyright (C) 2015, 2022 Rafael Vazquez
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

function errl2 = sp_l2_error (space, msh, u, uex)

  if (space.npatch ~= msh.npatch)
    error ('The number of patches does not coincide') 
  end

  error_l2 = zeros (msh.npatch, 1);
  for iptc = 1:msh.npatch
    if (space.ndof == numel(u))
  %     if (isempty (space.dofs_ornt))
      [Cpatch, Cpatch_cols] = sp_compute_Cpatch (space, iptc);
      u_ptc = Cpatch * u(Cpatch_cols);
  %     else
  %       u_ptc = u(space.gnum{iptc}) .* space.dofs_ornt{iptc}.';
  %     end
      error_l2(iptc) = sp_l2_error (space.sp_patch{iptc}, msh.msh_patch{iptc}, u_ptc, uex);
    else
      [Cpatch, Cpatch_cols] = sp_compute_Cpatch_vector (space, iptc, msh.rdim);
      u_ptc = Cpatch * u(Cpatch_cols);
      sp_vec = sp_vector (repmat (space.sp_patch(iptc), msh.rdim, 1), msh.msh_patch{iptc});
      error_l2(iptc) = sp_l2_error (sp_vec, msh.msh_patch{iptc}, u_ptc, uex);
    end
  end
  errl2 = sqrt (sum (error_l2 .* error_l2));

end
