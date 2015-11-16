% SP_L2_ERROR: Evaluate the error in L^2 norm.
%
%   [errl2, errl2_elem] = sp_l2_error (space, msh, u, uex);
%
% INPUT:
%
%     space:    structure representing the space of discrete functions (see sp_scalar/sp_evaluate_col)
%     msh:      structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%     u:        vector of dof weights
%     uex:      function handle to evaluate the exact solution
%
% OUTPUT:
%
%     errl2:  error in L^2 norm
%     errl2_elem:  error in L^2 norm, for each single element
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

function [errl2, errl2_elem] = sp_l2_error (sp, msh, u, uex)
  
  valu = sp_eval_msh (u, sp, msh);
  valu = reshape (valu, sp.ncomp, msh.nqn, msh.nel);

  for idir = 1:msh.rdim
    x{idir} = reshape (msh.geo_map(idir,:,:), msh.nqn*msh.nel, 1);
  end
  valex  = reshape (feval (uex, x{:}), sp.ncomp, msh.nqn, msh.nel);

  w = msh.quad_weights .* msh.jacdet;

  errl2_elem = sum (reshape (sum ((valu - valex).^2, 1), [msh.nqn, msh.nel]) .* w);

  errl2  = sqrt (sum (errl2_elem));
  errl2_elem  = sqrt (errl2_elem);
  
end
