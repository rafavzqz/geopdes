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

function [err_stress, err_stress_elem] = sp_stress_error (sp, msh, u, lambda_lame, mu_lame, graduex)

  stress_valu = sp_eval_msh (u, sp, msh, 'stress', lambda_lame, mu_lame);
  stress_valu = reshape (stress_valu, sp.ncomp, msh.rdim, msh.nqn, msh.nel);

  for idir = 1:msh.rdim
    x{idir} = reshape (msh.geo_map(idir,:,:), msh.nqn*msh.nel, 1);
  end
  grad_valex  = reshape (feval (graduex, x{:}), sp.ncomp, msh.rdim, msh.nqn, msh.nel);
  gradt_valex = permute (grad_valex, [2 1 3 4]);
  div_valex = 0;
  for ii = 1:msh.rdim
    div_valex = div_valex + grad_valex(ii,ii,:,:);
  end

  coeff_mu = reshape (mu_lame(x{:}), [1, 1, msh.nqn, msh.nel]);
  coeff_lambda = reshape (lambda_lame(x{:}), [1, 1, msh.nqn, msh.nel]);
  stress1 = bsxfun (@times, coeff_mu, (grad_valex + gradt_valex));
  div_valex = coeff_lambda .* div_valex;
  for ii = 1:msh.rdim
    stress2(ii,ii,:,:) = div_valex;
  end
  stress_valex = stress1 + stress2;

  w = msh.quad_weights .* msh.jacdet;

  err_stress_elem = sum (reshape (sum (sum ((stress_valu - stress_valex).^2, 1), 2), [msh.nqn, msh.nel]) .* w);
  err_stress = sqrt (sum (err_stress_elem));

  err_stress_elem = sqrt (err_stress_elem);
  
end
