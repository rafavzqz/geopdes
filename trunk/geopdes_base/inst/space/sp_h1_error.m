% SP_H1_ERROR: Evaluate the error in H^1 norm.
%
%   [toterr, errl2] = sp_h1_error (space, msh, u, uex, graduex);
%
% INPUT:
%
%     space:    structure representing the space of discrete functions (see sp_bspline_3d_phys)
%     msh:      structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%     u:        vector of dof weights
%     uex:      function handle to evaluate the exact solution
%     graduex:  function handle to evaluate the gradient of the exact solution
%
% OUTPUT:
%
%     toterr: error in H^1 norm
%     errl2:  error in L^2 norm
%
% Copyright (C) 2010 Carlo de Falco
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

function [toterr, errl2] = sp_h1_error (sp, msh, u, uex, graduex);

  ndir = size (msh.geo_map, 1);

  gradients = reshape (sp.shape_function_gradients, sp.ncomp, ndir, msh.nqn, sp.nsh_max, msh.nel);
  grad_num  = zeros (sp.ncomp, ndir, msh.nqn, msh.nel);
  for idir = 1:ndir
    for ish=1:sp.nsh_max
      for icmp=1:sp.ncomp
        grad_num(icmp, idir, :, :) = reshape (grad_num(icmp, idir, :, :), [msh.nqn, msh.nel]) + ...
          reshape (gradients (icmp, idir, :, ish, :), [msh.nqn, msh.nel]) .* ...
          repmat (u(sp.connectivity(ish,:)), 1, msh.nqn).';
      end
    end
  end

  
  for idir = 1:ndir
    x{idir} = reshape (msh.geo_map(idir,:,:), msh.nqn*msh.nel, 1);
  end
  grad_valex  = reshape (feval (graduex, x{:}), sp.ncomp, ndir, msh.nqn, msh.nel);

  w = msh.quad_weights .* msh.jacdet;
  errh1s = sum (sum (reshape (sum (sum ((grad_num - grad_valex).^2, 1), 2), [msh.nqn, msh.nel]) .* w));
  
  errl2  = sp_l2_error (sp, msh, u, uex);
  toterr = sqrt (errl2^2 + errh1s);
  
end

%!test
%! geometry = geo_load ({(@(PTS) PTS), (@(PTS) repmat (eye(2), [1, 1, size(PTS,2)]))});
%! knots = {[0 0 0 .5 1 1 1], [0 0 0 .33 .66 1 1 1]};
%! [qn, qw] = msh_set_quad_nodes (knots, msh_gauss_nodes ([3 3]));
%! msh = msh_2d_tensor_product (knots, qn, qw); 
%! msh = msh_push_forward_2d (msh, geometry);
%! sp  = sp_bspline_2d_phys (knots, [2 2], msh); 
%! u   = zeros (sp.ndof, 1);
%! uex     = @(x, y) zeros (size (x));
%! graduex = @(x, y) cat (1, zeros (size (x)), zeros (size (x)));
%! assert (sp_h1_error (sp, msh, u, uex, graduex), 0, eps);
