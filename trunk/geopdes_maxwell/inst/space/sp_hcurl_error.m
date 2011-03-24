% SP_HCURL_ERROR: Evaluate the error in H(curl) norm.
%
%   [toterr, errl2] = sp_hcurl_error (space, msh, u, uex, curluex);
%
% INPUT:
%
%     space:    structure representing the space of discrete functions (see sp_bspline_hcurl_3d)
%     msh:      structure containing the domain partition and the quadrature rule (see msh_push_forward_3d)
%     u:        vector of dof weights
%     uex:      function handle to evaluate the exact solution
%     curluex:  function handle to evaluate the curl of the exact solution
%
% OUTPUT:
%
%     toterr: error in H(curl) norm
%     errl2:  error in L^2 norm
%
% Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
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

function [toterr, errl2] = sp_hcurl_error (sp, msh, u, uex, curluex);

  w = msh.quad_weights(:) .* msh.jacdet(:);
  ndir = size (msh.geo_map, 1);
  for idir = 1:ndir
    x{idir} = reshape (msh.geo_map(idir,:,:), msh.nqn*msh.nel, 1);
  end

  curl_valex = feval (curluex, x{:});
  err_curl = 0;
  switch (ndir)
   case{2}
     valnum = zeros (msh.nqn, msh.nel);
     for ish = 1:sp.nsh_max
       valnum = valnum + ...
         reshape (sp.shape_function_curls(:, ish, :), msh.nqn, msh.nel) .* ...
         u(repmat(sp.connectivity(ish, :), msh.nqn, 1));
     end
     err_curl = err_curl + sum((valnum(:) - curl_valex).^2 .* w);
   case{3}
    for idir = 1:ndir
      valnum = zeros (msh.nqn, msh.nel);
      for ish=1:sp.nsh_max
        valnum = valnum + ...
          reshape (sp.shape_function_curls(idir, :, ish, :), msh.nqn, msh.nel) .* ...
	  reshape (u(repmat(sp.connectivity(ish,:), msh.nqn, 1)), msh.nqn, msh.nel);
      end
      valex = curl_valex(idir, :);
      err_curl = err_curl + sum ((valnum(:) - valex(:)).^2 .* w);
    end
  end

  errl2  = sp_l2_error (sp, msh, u, uex);
  toterr = sqrt (errl2^2 + err_curl);

end
