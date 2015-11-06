% SP_H1_ERROR: Evaluate the error in H^1 norm.
%
%   [toterr, errl2] = sp_h1_error (space, msh, u, uex, graduex)
%
% INPUT:
%
%    space:   object defining the space of discrete functions (see sp_bspline_2d)
%    msh:     object defining the domain partition and the quadrature rule (see msh_2d)
%    u:       vector of dof weights
%    uex:     function handle to evaluate the exact solution
%    graduex: function handle to evaluate the gradient of the exact solution
%
% OUTPUT:
%
%     toterr: error in H^1 norm
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

function [errh1, errl2] = sp_h1_error (sp, msh, u, uex, graduex)

  ndim = numel (msh.qn);

  errl2 = 0;
  errh1s = 0;

  valu = zeros (sp.ncomp, msh.nqn, msh.nelcol);
  grad_valu = zeros (sp.ncomp, ndim, msh.nqn, msh.nelcol);
  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (sp, msh_col, 'gradient', true);

    uc_iel = zeros (size (sp_col.connectivity));
    uc_iel(sp_col.connectivity~=0) = ...
          u(sp_col.connectivity(sp_col.connectivity~=0));

    weight = repmat (reshape (uc_iel, [1, sp_col.nsh_max, msh_col.nel]), ...
                                  [msh.nqn, 1, 1]);

    sp_col.shape_functions = reshape (sp_col.shape_functions, ...
           sp.ncomp, msh_col.nqn, sp_col.nsh_max, msh_col.nel);
    sp_col.shape_function_gradients = reshape (sp_col.shape_function_gradients, ...
           sp.ncomp, ndim, msh_col.nqn, sp_col.nsh_max, msh_col.nel);

    for icmp = 1:sp.ncomp
      for idim = 1:ndim
        grad_valu(icmp,idim,:,:) = sum (weight .* reshape (sp_col.shape_function_gradients(icmp,idim,:,:,:), ...
                            msh_col.nqn, sp_col.nsh_max, msh_col.nel), 2);
      end
      valu(icmp,:,:) = sum (weight .* reshape (sp_col.shape_functions(icmp,:,:,:), ...
                            msh_col.nqn, sp_col.nsh_max, msh_col.nel), 2);
    end

    for idim = 1:ndim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end
    w = msh_col.quad_weights .* msh_col.jacdet;

    valex  = reshape (feval (uex, x{:}), sp.ncomp, msh_col.nqn, msh_col.nel);
    grad_valex  = reshape (feval (graduex, x{:}), sp.ncomp, ndim, msh_col.nqn, msh_col.nel);

    errh1s = errh1s + sum (sum (reshape (sum (sum ...
       ((grad_valu - grad_valex).^2, 1), 2), [msh_col.nqn, msh_col.nel]) .* w));

    erraux = sum ((valu - valex).^2, 1);
    errl2  = errl2 + sum (w(:) .* erraux(:));

  end

  errh1 = sqrt (errl2 + errh1s);
  errl2 = sqrt (errl2);
  
end

