% SP_H1_ERROR: Evaluate the error in H^1 norm.
%
%   [errh1, errl2, errh1s, errh1_elem, errl2_elem, errh1s_elem] = sp_h1_error (space, msh, u, uex, graduex);
%
% INPUT:
%
%     space:    structure representing the space of discrete functions (see sp_scalar/sp_evaluate_col)
%     msh:      structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%     u:        vector of dof weights
%     uex:      function handle to evaluate the exact solution
%     graduex:  function handle to evaluate the gradient of the exact solution
%
% OUTPUT:
%
%     errh1:  error in H^1 norm
%     errl2:  error in L^2 norm
%     errh1s: error in H^1 seminorm
%     errh1_elem:  error in H^1 norm, for each single element
%     errl2_elem:  error in L^2 norm, for each single element
%     errh1s_elem: error in H^1 seminorm, for each single element
%
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2011 Rafael Vazquez
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
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

function [errh1, errl2, errh1s, errh1_elem, errl2_elem, errh1s_elem] = sp_h1_error (sp, msh, u, uex, graduex)

  grad_valu = sp_eval_msh (u, sp, msh, 'gradient');
  grad_valu = reshape (grad_valu, sp.ncomp, msh.rdim, msh.nqn, msh.nel);

  for idir = 1:msh.rdim
    x{idir} = reshape (msh.geo_map(idir,:,:), msh.nqn*msh.nel, 1);
  end
  grad_valex  = reshape (feval (graduex, x{:}), sp.ncomp, msh.rdim, msh.nqn, msh.nel);

  w = msh.quad_weights .* msh.jacdet;

  [errl2, errl2_elem] = sp_l2_error (sp, msh, u, uex);
  errh1s_elem = sum (reshape (sum (sum ((grad_valu - grad_valex).^2, 1), 2), [msh.nqn, msh.nel]) .* w);
  errh1s = sqrt (sum (errh1s_elem));

  errh1  = sqrt (errl2^2 + errh1s^2);

  errh1_elem  = sqrt (errl2_elem.^2 + errh1s_elem);
  errh1s_elem = sqrt (errh1s_elem);
  
end
