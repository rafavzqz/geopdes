% OP_SU_EV_MP: assemble the matrix A = [a(i,j)], a(i,j) = 1/2 (sigma (u_j), epsilon (v_i)), in a multipatch domain.
%
%   mat = op_su_ev_mp (spu, spv, msh, lambda, mu, [patches]);
%
% INPUT:
%    
%   spu:     object representing the space of trial functions (see sp_multipatch)
%   spv:     object representing the space of test functions (see sp_multipatch)
%   msh:     object that defines the domain partition and the quadrature rule (see msh_multipatch)
%   lambda, mu: function handles to compute the Lame' coefficients
%   patches: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%   mat:    assembled matrix
% 
% Copyright (C) 2015, 2022 Rafael Vazquez
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

    [Cpatch_u, Cpatch_cols_u] = sp_compute_Cpatch (spu, iptc);
    [Cpatch_v, Cpatch_cols_v] = sp_compute_Cpatch (spv, iptc);
    
    Cpatch_u = repmat ({Cpatch_u}, 1, msh.rdim);
    Cpatch_v = repmat ({Cpatch_v}, 1, msh.rdim);
    Cpatch_vector_spu = blkdiag (Cpatch_u{:});
    Cpatch_vector_spv = blkdiag (Cpatch_v{:});
    
    cols_u = []; cols_v = [];
    for icomp = 1:msh.rdim
      cols_u = union (cols_u, (icomp-1)*spu.ndof + Cpatch_cols_u);
      cols_v = union (cols_v, (icomp-1)*spv.ndof + Cpatch_cols_v);
    end
    
    A(cols_v,cols_u) = ...
      A(cols_v,cols_u) + Cpatch_vector_spv.' * Ap * Cpatch_vector_spu;
  end

end