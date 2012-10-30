% TSPLINE_H1_ERROR: Evaluate the error in H^1 norm.
%
%   [errh1, errl2] = tspline_h1_error (space, msh, u, uex, graduex);
%
% INPUT:
%
%     space:    structure representing the space of discrete functions (see tspline_mesh_space)
%     msh:      structure containing the domain partition and the quadrature rule (see tspline_mesh_space)
%     u:        vector of dof weights
%     uex:      function handle to evaluate the exact solution
%     graduex:  function handle to evaluate the gradient of the exact solution
%
% OUTPUT:
%
%     errh1: error in H^1 norm
%     errl2:  error in L^2 norm
%
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2012 Rafael Vazquez
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

function [errh1, errl2] = tspline_h1_error (sp, msh, u, uex, graduex)

  rdim = size (msh.geo_map_jac, 1);
  ndim = size (msh.geo_map_jac, 2);
  
  valu = tspline_sp_eval_msh (u, sp, msh, 'value');
  grad_valu = tspline_sp_eval_msh (u, sp, msh, 'gradient');

  valu = reshape (valu, sp.ncomp, msh.nqn, msh.nel);
  grad_valu = reshape (grad_valu, sp.ncomp, ndim, msh.nqn, msh.nel);

% Compute the exact solution in the quadrature points
  for idir = 1:rdim
    x{idir} = reshape (msh.geo_map(idir,:,:), msh.nqn*msh.nel, 1);
  end
  valex = reshape (feval (uex, x{:}), sp.ncomp, msh.nqn, msh.nel);
  grad_valex  = reshape (feval (graduex, x{:}), sp.ncomp, ndim, msh.nqn, msh.nel);

% Evaluate the integrals to compute the error
  w = msh.quad_weights .* msh.jacdet;
  errh1s = sum (sum (reshape (sum (sum ((grad_valu - grad_valex).^2, 1), 2), [msh.nqn, msh.nel]) .* w));
  errl2 = sum ((valu - valex).^2, 1);
  errl2 = sqrt (sum (errl2(:) .* w(:)));

  errh1 = sqrt (errl2^2 + errh1s);
  
end
