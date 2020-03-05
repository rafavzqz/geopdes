% SP_HDIV_ERROR: Evaluate the error in H(div) norm.
%
%   [errhdiv, errl2, errdiv] = sp_hdiv_error (space, msh, u, uex, divuex)
%
% INPUT:
%
%    space:   object defining the space of discrete functions (see sp_vector)
%    msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%    u:       vector of dof weights
%    uex:     function handle to evaluate the exact solution
%    divuex:  function handle to evaluate the divergence of the exact solution
%
% OUTPUT:
%
%     errhdiv: error in H(div) norm
%     errl2:   error in L^2 norm
%     errdiv:  error of the divergence in L^2 norm
%
% Copyright (C) 2020 Rafael Vazquez
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

function [errhdiv, errl2, errdiv] = sp_hdiv_error (space, msh, u, uex, divuex)

  if (numel(u) ~= space.ndof)
    error ('Wrong size of the vector of degrees of freedom')
  end

  errl2 = 0; errdiv = 0;
  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col, 'value', true, 'divergence', true);
    
    [~, err_l2, err_div] = sp_hdiv_error (sp_col, msh_col, u, uex, divuex);
    
    errdiv = errdiv + err_div.^2;
    errl2 = errl2 + err_l2.^2;
  end
  
  errhdiv = sqrt (errl2 + errdiv);
  errl2 = sqrt (errl2);
  errdiv = sqrt (errdiv);

end
