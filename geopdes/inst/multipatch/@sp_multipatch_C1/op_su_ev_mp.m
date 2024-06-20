% OP_SU_EV_MP: assemble the matrix A = [a(i,j)], a(i,j) = 1/2 (sigma (u_j), epsilon (v_i)), in a multipatch domain.
%
%   mat = op_su_ev_mp (spu, spv, msh, lambda, mu, [patches]);
%
% INPUT:
%    
%   spu:     object representing the space of trial functions (see sp_multipatch_C1)
%   spv:     object representing the space of test functions (see sp_multipatch_C1)
%   msh:     object that defines the domain partition and the quadrature rule (see msh_multipatch)
%   lambda, mu: function handles to compute the Lame' coefficients
%   patches: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%   mat:    assembled matrix
% 
% Copyright (C) 2015, 2022, 2023 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function A = op_su_ev_mp (spu, spv, msh, lambda, mu, patch_list)

  if (nargin < 6)
    patch_list = 1:msh.npatch;
  end

  if ((spu.npatch ~= spv.npatch) || (spu.npatch ~= msh.npatch))
    error ('op_su_ev_mp: the number of patches does not coincide')
  end
  
  A = sparse (msh.rdim*spv.ndof, msh.rdim*spu.ndof);
  
  for iptc = patch_list
    Ap = op_su_ev_tp (spu.sp_patch{iptc}, spv.sp_patch{iptc}, msh.msh_patch{iptc}, lambda, mu);

    [Cpatch_u, cols_u] = sp_compute_Cpatch_vector (spu, iptc, msh.rdim);
    [Cpatch_v, cols_v] = sp_compute_Cpatch_vector (spv, iptc, msh.rdim);
        
    A(cols_v,cols_u) = A(cols_v,cols_u) + Cpatch_v.' * Ap * Cpatch_u;
  end

end