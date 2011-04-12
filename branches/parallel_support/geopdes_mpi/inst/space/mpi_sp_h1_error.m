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
% Copyright (C) 2011 Andrea Bressan
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

function [toterr, errl2] = mpi_sp_h1_error (comm, sp, msh, u, uex, graduex);
  [loc_h1, loc_l2] = sp_h1_error (sp, msh, u, uex, graduex);
  errl2s           = MPI_Allreduce(loc_l2^2, 'mpi_sum', comm);
  errh1s           = MPI_Allreduce(loc_h1^2, 'mpi_sum', comm);
  errl2  = sqrt (errl2s);
  toterr = sqrt (errh1s);
end
