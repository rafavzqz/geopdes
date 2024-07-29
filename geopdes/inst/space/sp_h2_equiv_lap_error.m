% SP_H2_ERROR: Evaluate the error in H^2 equivalent seminorm with laplacian, H^1 and L^2 norms.
%
%   [errh2, errh1, errl2, errh2s, errh1s, errh2_elem, errh1_elem, errl2_elem, errh2s_elem, errh1s_elem] = sp_h2_error (sp, msh, u, uex, graduex, lapuex)
%
% INPUT:
%
%    space:   structure representing the space of discrete functions (see sp_scalar/sp_evaluate_col)
%    msh:     structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%    u:       vector of dof weights
%    uex:     function handle to evaluate the exact solution
%    graduex: function handle to evaluate the gradient of the exact solution
%    lapuex:  function handle to evaluate the laplacian of the exact solution
%
% OUTPUT:
%
%     errh2:  error in (equivalent) H^2 norm
%     errh1:  error in H^1 norm
%     errl2:  error in L^2 norm
%     errh2s: error in H^2 seminorm (involving only second dervatives)
%     errh1s: error in H^1 seminorm
%     errh2_elem:  error in H^2 norm, for each single element
%     errh1_elem:  error in H^1 norm, for each single element
%     errl2_elem:  error in L^2 norm, for each single element
%     errh2s_elem: error in H^2 seminorm, for each single element
%     errh1s_elem: error in H^1 seminorm, for each single element
%
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2011, 2016 Rafael Vazquez
% Copyright (C) 2015, 2016 Viacheslav Balobanov
% Copyright (C) 2022 Andrea Farahat
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

function [errh2, errh1, errl2, errh2s, errh1s, errh2_elem, errh1_elem, errl2_elem, errh2s_elem, errh1s_elem] = ...
    sp_h2_equiv_lap_error (sp, msh, u, uex, graduex, lapuex)
  
  lap_valu = sp_eval_msh (u, sp, msh, 'laplacian');
  lap_valu = reshape (lap_valu, sp.ncomp, msh.nqn, msh.nel);

  for idir = 1:msh.rdim
    x{idir} = reshape (msh.geo_map(idir,:,:), msh.nqn*msh.nel, 1);
  end
  lap_valex  = reshape (feval (lapuex, x{:}), sp.ncomp, msh.nqn, msh.nel);

  w = msh.quad_weights .* msh.jacdet;

  [errh1, errl2, errh1s, errh1_elem, errl2_elem, errh1s_elem] = sp_h1_error (sp, msh, u, uex, graduex);

  errh2s_elem = sum (reshape (sum ((lap_valu - lap_valex).^2, 1), [msh.nqn, msh.nel]) .*w);
  errh2s = sqrt (sum (errh2s_elem));

  errh2 = sqrt (errh1^2 + errh2s^2);
  errh2_elem = sqrt (errh1_elem.^2 + errh2s_elem);

end