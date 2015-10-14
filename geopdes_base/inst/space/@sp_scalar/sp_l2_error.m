% SP_L2_ERROR: Evaluate the error in L^2 norm.
%
%   errl2 = sp_l2_error (space, msh, u, uex)
%
% INPUT:
%
%   space: object defining the space of discrete functions (see sp_scalar)
%   msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%   u:     vector of dof weights
%   uex:   function handle to evaluate the exact solution
%
% OUTPUT:
%
%     errl2:  error in L^2 norm
%
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2011, 2015 Rafael Vazquez
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

  errl2 = 0;
  
  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col, 'value', true, 'gradient', false);
    
    errl2 = errl2 + (sp_l2_error (sp_col, msh_col, u, uex)).^2;
  end
  
  errl2 = sqrt (errl2);

end
