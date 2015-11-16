% SP_HCURL_ERROR: Evaluate the error in H(curl) norm.
%
%   [errhcurl, errl2, errcurl] = sp_hcurl_error (space, msh, u, uex, curluex);
%
% INPUT:
%
%    space:   object defining the space of discrete functions (see sp_multipatch)
%    msh:     object defining the domain partition and the quadrature rule (see msh_multipatch)
%    u:       vector of dof weights
%    uex:     function handle to evaluate the exact solution
%    curluex: function handle to evaluate the curl of the exact solution
%
% OUTPUT:
%
%     errhcurl: error in H(curl) norm
%     errl2:    error in L^2 norm
%     errcurl:  error of the curl in L^2 norm
%
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

function [errhcurl, errl2, errcurl] = sp_hcurl_error (space, msh, u, uex, curluex)

  if (space.npatch ~= msh.npatch)
    error ('The number of patches does not coincide') 
  end

  for iptc = 1:msh.npatch
    if (isempty (space.dofs_ornt))
      u_ptc = u(space.gnum{iptc});
    else
      u_ptc = u(space.gnum{iptc}) .* space.dofs_ornt{iptc}.';
    end
    [error_hcurl(iptc), error_l2(iptc), error_curl(iptc)] = ...
        sp_hcurl_error (space.sp_patch{iptc}, msh.msh_patch{iptc}, u_ptc, uex, curluex);
  end
  errl2 = sqrt (sum (error_l2 .* error_l2));
  errhcurl = sqrt (sum (error_hcurl .* error_hcurl));
  errcurl = sqrt (sum (error_curl .* error_curl));

end
