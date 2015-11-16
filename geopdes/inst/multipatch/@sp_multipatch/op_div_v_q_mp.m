% OP_DIV_V_Q_MP: assemble the matrix B = [b(i,j)], b(i,j) = (q_i, div v_j), in a multipatch domain.
%
%   mat = op_div_v_q_mp (spv, spq, msh, [patches])
%
% INPUT: 
%
%   spv:     object representing the space of trial functions for the velocity (see sp_multipatch))
%   spq:     object representing the space of test functions for the pressure (see sp_multipatch)
%   msh:     object that defines the domain partition and the quadrature rule (see msh_multipatch)
%   patches: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT: 
%
%   mat:    assembled matrix 
% 
% Copyright (C) 2015 Rafael Vazquez
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

function A = op_div_v_q_mp (spv, spq, msh, patch_list)

  if (nargin < 4)
    patch_list = 1:msh.npatch;
  end

  if ((spv.npatch ~= spq.npatch) || (spv.npatch ~= msh.npatch))
    error ('op_div_v_q_mp: the number of patches does not coincide')
  end
  
  ncounter = 0;
  for iptc = patch_list
    [rs, cs, vs] = op_div_v_q_tp (spv.sp_patch{iptc}, spq.sp_patch{iptc}, msh.msh_patch{iptc});
    rows(ncounter+(1:numel (rs))) = spq.gnum{iptc}(rs);
    cols(ncounter+(1:numel (rs))) = spv.gnum{iptc}(cs);

    if (~isempty (spv.dofs_ornt))
      vs = vs .* spv.dofs_ornt{iptc}(cs)';
    end
    
    vals(ncounter+(1:numel (rs))) = vs;
    ncounter = ncounter + numel (rs);
  end

  A = sparse (rows, cols, vals, spq.ndof, spv.ndof);
  clear rows cols vals rs cs vs

end