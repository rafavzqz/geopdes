% SP_HDIV_ERROR: Evaluate the error in H(div) norm.
%
%   [errhdiv, errl2, errhdivs, errhdiv_elem, errl2_elem, errhdivs_elem] = sp_hdiv_error (space, msh, u, uex, divuex);
%
% INPUT:
%
%    space:   object defining the space of discrete functions (see sp_vector)
%    msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%    u:       vector of dof weights
%    uex:     function handle to evaluate the exact solution
%    divuex:  function handle to evaluate the div of the exact solution
%
% OUTPUT:
%
%     errhdiv:       error in H(div) norm
%     errl2:         error in L^2 norm
%     errhdivs:      error of the div in L^2 norm
%     errhdiv_elem:  error in H(div) norm, for each single element
%     errl2_elem:    error in L^2 norm, for each single element
%     errhdivs_elem: error of the div in L^2 norm, for each single element
%
% Copyright (C) 2020 Riccardo Puppi
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
% <http://www.gnu.org/licenses/>.function [errhdiv, errl2, errhdivs, errhdiv_elem, errl2_elem, errhdivs_elem] = sp_hdiv_error (sp, msh, u, uex, divuex)

function [errhdiv, errl2, errhdivs, errhdiv_elem, errl2_elem, errhdivs_elem] = sp_hdiv_error (sp, msh, u, uex, divuex)

  div_valu = sp_eval_msh (u, sp, msh, 'divergence');
  div_valu = reshape (div_valu, msh.nqn, msh.nel);

  for idir = 1:msh.rdim
    x{idir} = reshape (msh.geo_map(idir,:,:), msh.nqn*msh.nel, 1);
  end
  div_valex  = reshape (feval (divuex, x{:}), msh.nqn, msh.nel);

  w = msh.quad_weights .* msh.jacdet;

  [errl2, errl2_elem] = sp_l2_error (sp, msh, u, uex);
  errhdivs_elem = sum ((div_valu - div_valex).^2 .* w);
  errhdivs = sqrt (sum (errhdivs_elem));

  errhdiv  = sqrt (errl2^2 + errhdivs^2);

  errhdiv_elem  = sqrt (errl2_elem.^2 + errhdivs_elem);
  errhdivs_elem = sqrt (errhdivs_elem);
  
end
