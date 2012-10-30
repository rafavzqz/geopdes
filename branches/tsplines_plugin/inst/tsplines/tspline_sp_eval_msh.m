% TSPLINE_SP_EVAL_MSH: Evaluate a function, given by its degrees of freedom, at the points given by a msh object.
%
%   [eu, F] = tspline_sp_eval_msh (u, space, msh, [option]);
%
% INPUT:
%
%     u:         vector of dof weights
%     space:     structure defining the discrete space (see tspline_space_mesh)
%     msh:       structure defining the points where to evaluate (see tspline_space_mesh)
%     option:    accepted options are 'value' (default) and 'gradient'.
%
% OUTPUT:
%
%     eu: the function evaluated in the given points 
%     F:  grid points in the physical domain
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

function [eu, F] = tspline_sp_eval_msh (u, sp, msh, option)

  if (nargin < 4)
    option = 'value';
  end

  ndim = size (msh.geo_map_jac, 2);

% Write the values of u in an array with the same size of connectivity
  u_conn = zeros (size (sp.connectivity));
  u_conn(sp.connectivity~=0) = u(sp.connectivity(sp.connectivity~=0));

  weight = reshape (u_conn, [1, sp.nsh_max, msh.nel]);
  clear u_conn


  if (strcmpi (option, 'value'))
    eu = zeros (sp.ncomp, msh.nqn, msh.nel);
% Reshape the function values, to use the same code for scalars and vectors
    sp.shape_functions = reshape (sp.shape_functions, ...
                                  sp.ncomp, msh.nqn, sp.nsh_max, msh.nel);

    for icmp = 1:sp.ncomp
      eu(icmp,:,:) = sum (bsxfun (@times, weight, ...
          reshape (sp.shape_functions(icmp,:,:,:), msh.nqn, sp.nsh_max, msh.nel)), 2);
    end

  elseif (strcmpi (option, 'gradient'))
    eu = zeros (sp.ncomp, ndim, msh.nqn, msh.nel);
% Reshape the function values, to use the same code for scalars and vectors
    sp.shape_function_gradients = reshape (sp.shape_function_gradients, ...
                                  sp.ncomp, ndim, msh.nqn, sp.nsh_max, msh.nel);

    for icmp = 1:sp.ncomp
      for idim = 1:ndim
        eu(icmp,idim,:,:) = sum (bsxfun (@times, weight, ...
                      reshape (sp.shape_function_gradients(icmp,idim,:,:,:), ...
                      msh.nqn, sp.nsh_max, msh.nel)), 2);
      end
    end
  end

  if (sp.ncomp == 1)
    saux = size (eu);
    eu = reshape (eu, saux(2:end));
  end
  F = msh.geo_map;
  
end
