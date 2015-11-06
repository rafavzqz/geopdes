% SP_HCURL_ERROR: Evaluate the error in H(curl) norm.
%
%   [toterr, errl2] = sp_hcurl_error (space, msh, u, uex, curluex);
%
% INPUT:
%
%    space:   object defining the space of discrete functions (see sp_vector_2d_curl_transform)
%    msh:     object defining the domain partition and the quadrature rule (see msh_2d)
%    u:       vector of dof weights
%    uex:     function handle to evaluate the exact solution
%    curluex: function handle to evaluate the curl of the exact solution
%
% OUTPUT:
%
%     toterr: error in H(curl) norm
%     errl2:  error in L^2 norm
%
% Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
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

function [toterr, errl2] = sp_hcurl_error (sp, msh, u, uex, curluex);

  ndim = numel (msh.qn);

  errl2 = 0;
  err_curl = 0;

  switch (ndim)
    case {2}
      valu = zeros (sp.ncomp, msh.nqn, msh.nelcol);

      for iel = 1:msh.nel_dir(1)
        msh_col = msh_evaluate_col (msh, iel);
        sp_col  = sp_evaluate_col (sp, msh_col, 'curl', true);

        uc_iel = zeros (size (sp_col.connectivity));
        uc_iel(sp_col.connectivity~=0) = ...
          u(sp_col.connectivity(sp_col.connectivity~=0));

        weight = repmat (reshape (uc_iel, [1, sp_col.nsh_max, msh_col.nel]), ...
                                  [msh_col.nqn, 1, 1]);
        sp_col.shape_functions = reshape (sp_col.shape_functions, ...
               sp.ncomp, msh_col.nqn, sp_col.nsh_max, msh_col.nel);
        sp_col.shape_function_curls = reshape (sp_col.shape_function_curls, ...
           msh_col.nqn, sp_col.nsh_max, msh_col.nel);

        for icmp = 1:sp.ncomp
          valu(icmp,:,:) = sum (weight .* reshape (sp_col.shape_functions(icmp,:,:), ...
                            msh_col.nqn, sp_col.nsh_max, msh_col.nel), 2);
        end
        curl_valu = sum (weight .* sp_col.shape_function_curls, 2);
        curl_valu = reshape (curl_valu, msh_col.nqn, msh_col.nel);

        for idim = 1:ndim
          x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
        end
        w = msh_col.quad_weights .* msh_col.jacdet;

        valex  = reshape (feval (uex, x{:}), sp.ncomp, msh_col.nqn, msh_col.nel);
        curl_valex  = reshape (feval (curluex, x{:}), msh_col.nqn, msh_col.nel);

        err_curl = err_curl + sum (w(:) .* (curl_valu(:) - curl_valex(:)).^2);

        erraux = sum ((valu - valex).^2, 1);
        errl2  = errl2 + sum (w(:) .* erraux(:));
      end

    case {3}
      valu = zeros (sp.ncomp, msh.nqn, msh.nelcol);
      curl_valu = zeros (sp.ncomp, msh.nqn, msh.nelcol);

      for iel = 1:msh.nel_dir(1)
        msh_col = msh_evaluate_col (msh, iel);
        sp_col  = sp_evaluate_col (sp, msh_col, 'curl', true);

        uc_iel = zeros (size (sp_col.connectivity));
        uc_iel(sp_col.connectivity~=0) = ...
          u(sp_col.connectivity(sp_col.connectivity~=0));

        weight = repmat (reshape (uc_iel, [1, sp_col.nsh_max, msh_col.nel]), ...
                                  [msh_col.nqn, 1, 1]);
        sp_col.shape_functions = reshape (sp_col.shape_functions, ...
               sp.ncomp, msh_col.nqn, sp_col.nsh_max, msh_col.nel);
        sp_col.shape_function_curls = reshape (sp_col.shape_function_curls, ...
           sp.ncomp, msh_col.nqn, sp_col.nsh_max, msh_col.nel);

        for icmp = 1:sp.ncomp
          curl_valu(icmp,:,:) = sum (weight .* reshape (sp_col.shape_function_curls(icmp,:,:), ...
                            msh_col.nqn, sp_col.nsh_max, msh_col.nel), 2);
          valu(icmp,:,:) = sum (weight .* reshape (sp_col.shape_functions(icmp,:,:), ...
                            msh_col.nqn, sp_col.nsh_max, msh_col.nel), 2);
        end

        for idim = 1:ndim
          x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
        end
        w = msh_col.quad_weights .* msh_col.jacdet;

        valex  = reshape (feval (uex, x{:}), sp.ncomp, msh_col.nqn, msh_col.nel);
        curl_valex  = reshape (feval (curluex, x{:}), sp.ncomp, msh_col.nqn, msh_col.nel);

        err_curl = err_curl + sum (sum (reshape (sum ...
         ((curl_valu - curl_valex).^2, 1), [msh_col.nqn, msh_col.nel]) .* w));

        erraux = sum ((valu - valex).^2, 1);
        errl2  = errl2 + sum (w(:) .* erraux(:));
      end

  end

  toterr = sqrt (errl2 + err_curl);
  errl2 = sqrt (errl2);
  
end
