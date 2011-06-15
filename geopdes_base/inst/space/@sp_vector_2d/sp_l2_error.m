% SP_L2_ERROR: Evaluate the error in L^2 norm.
%
%   errl2 = sp_l2_error (space, msh, u, uex);
%
% INPUT:
%
%    space:  class defining the space of discrete functions (see sp_vector_2d)
%    msh:    class defining the domain partition and the quadrature rule (see msh_2d)
%    u:      vector of dof weights
%    uex:    function handle to evaluate the exact solution
%
% OUTPUT:
%
%     errl2:  error in L^2 norm
%
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2011 Rafael Vazquez
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

function errl2 = sp_l2_error (sp, msh, u, uex)

  errl2 = 0;

  valu = zeros (2, msh.nqn, msh.nelv);
  for iel = 1:msh.nelu
    [sp_col, elem_list] = sp_evaluate_col (sp, msh, iel, 'gradient', false);

    uc_iel = zeros (size (sp_col.connectivity));
    uc_iel(sp_col.connectivity~=0) = ...
          u(sp_col.connectivity(sp_col.connectivity~=0));

    weight = repmat (reshape (uc_iel, [1, sp_col.nsh_max, msh.nelv]), ...
                                  [msh.nqn, 1, 1]);

    valu(1,:,:) = sum (weight .* reshape (sp_col.shape_functions(1,:,:,:), ...
                          msh.nqn, sp_col.nsh_max, msh.nelv), 2);
    valu(2,:,:) = sum (weight .* reshape (sp_col.shape_functions(2,:,:,:), ...
                          msh.nqn, sp_col.nsh_max, msh.nelv), 2);

    w = msh.quad_weights(:, elem_list) .* msh.jacdet(:, elem_list);
    x = msh.geo_map(1, :, elem_list);
    y = msh.geo_map(2, :, elem_list);

    valex = reshape (feval (uex, x(:), y(:)), sp.ncomp, [], msh.nelv);

    erraux = sum ((valu - valex).^2, 1);
    errl2  = errl2 + sum (w(:) .* erraux(:));
  end

  errl2 = sqrt (errl2);
  
end

