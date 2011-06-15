% SP_L2_ERROR: Evaluate the error in L^2 norm.
%
%   errl2 = sp_l2_error (space, msh, u, uex, graduex);
%
% INPUT:
%
%    space:    class defining the space of discrete functions (see sp_nurbs_3d)
%    msh:      class defining the domain partition and the quadrature rule (see msh_push_forward_3d)
%    u:        vector of dof weights
%    uex:      function handle to evaluate the exact solution
%    graduex:  function handle to evaluate the gradient of the exact solution
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

function errl2 = sp_l2_error (sp, msh, u, uex);

  errl2 = 0;

  nel_col = msh.nelv * msh.nelw;
  for iel = 1:msh.nelu
    [sp_col, elem_list] = sp_evaluate_col (sp, msh, iel, 'gradient', false);

    uc_iel = zeros (size (sp_col.connectivity));
    uc_iel(sp_col.connectivity~=0) = ...
          u(sp_col.connectivity(sp_col.connectivity~=0));

    weight = repmat (reshape (uc_iel, [1, sp_col.nsh_max, nel_col]), ...
                                  [msh.nqn, 1, 1]);

    valu = squeeze (sum (weight .* sp_col.shape_functions, 2));
 
    w = msh.quad_weights(:, elem_list) .* msh.jacdet(:, elem_list);
    x = msh.geo_map(1, :, elem_list);
    y = msh.geo_map(2, :, elem_list);
    z = msh.geo_map(3, :, elem_list);

    valex  = reshape (feval (uex, x, y, z), msh.nqn, nel_col);

    errl2  = errl2 + sum (w(:) .* (valu(:) - valex(:)).^2);

  end

  errl2 = sqrt (errl2);
  
end

