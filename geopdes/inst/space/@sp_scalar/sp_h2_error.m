% SP_H2_ERROR: Evaluate the error in H^2 norm, H^1 and L^2 norms.
%
%   [errh2, errh1, errl2, errh2s, errh1s] = sp_h2_error (sp, msh, u, uex, graduex, hessuex)
%
% INPUT:
%
%    space:   object defining the space of discrete functions (see sp_scalar)
%    msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%    u:       vector of dof weights
%    uex:     function handle to evaluate the exact solution
%    graduex: function handle to evaluate the gradient of the exact solution
%    hessuex: function handle to evaluate the hessian of the exact solution
%
% OUTPUT:
%
%     errh2: error in H^2 norm
%     errh1: error in H^1 norm
%     errl2: error in L^2 norm
%     errh2s: error in H^2 seminorm
%     errh1s: error in H^1 seminorm
%
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2011, 2016, 2017 Rafael Vazquez
% Copyright (C) 2015, 2016 Viacheslav Balobanov
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

function [errh2, errh1, errl2, errh2s, errh1s] = sp_h2_error (sp, msh, u, uex, graduex, hessuex)

  if (numel(u) ~= sp.ndof)
    error ('Wrong size of the vector of degrees of freedom')
  end

  errl2 = 0;
  errh1s = 0;
  errh2s = 0;

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (sp, msh_col, 'gradient', true, 'hessian', true);

    [~, ~, err_l2, err_h2s, err_h1s] = sp_h2_error (sp_col, msh_col, u, uex, graduex, hessuex);
    
    errh2s = errh2s + err_h2s.^2;
    errh1s = errh1s + err_h1s.^2;
    errl2 = errl2 + err_l2.^2;
  end
  
  errh2 = sqrt (errl2 + errh1s + errh2s);
  errh1 = sqrt (errl2 + errh1s);
  errl2 = sqrt (errl2);
  
  errh2s = sqrt (errh2s);
  errh1s = sqrt (errh1s);
  
end
