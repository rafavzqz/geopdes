% SP_HCURL_ERROR: Evaluate the error in H(curl) norm.
%
%   [errhcurl, errl2, errcurl] = sp_hcurl_error (space, msh, u, uex, curluex);
%
% INPUT:
%
%    space:   object defining the space of discrete functions (see sp_vector)
%    msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%    u:       vector of dof weights
%    uex:     function handle to evaluate the exact solution
%    curluex: function handle to evaluate the curl of the exact solution
%
% OUTPUT:
%
%     errhcurl: error in H(curl) norm
%     errl2:    error in L^2 norm
%     errcurl:  error of the curl in L^2 norm
%
% Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
% Copyright (C) 2011, 2015 Rafael Vazquez
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

function [errhcurl, errl2, errcurl] = sp_hcurl_error (space, msh, u, uex, curluex)

  errl2 = 0;
  errcurl = 0;

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col, 'value', true, 'curl', true);

    [~, err_l2, err_curl] = sp_hcurl_error (sp_col, msh_col, u, uex, curluex);

    errcurl = errcurl + err_curl.^2;
    errl2 = errl2 + err_l2.^2;
  end

  errhcurl = sqrt (errl2 + errcurl);
  errl2 = sqrt (errl2);
  errcurl = sqrt (errcurl);
  
end
