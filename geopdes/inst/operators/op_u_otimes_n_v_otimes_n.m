% OP_U_OTIMES_N_V_OTIMES_N: assemble the matrix A = [A(i,j)], A(i,j) = (epsilon (u \otimes n)_j, (v \otimes n)_i ).
%
%   A = op_u_otimes_n_v_otimes_n (spu, spv, msh, coeff);
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_vector/sp_eval_boundary_side)
%   spv:   structure representing the space of test functions (see sp_vector/sp_eval_boundary_side)
%   msh:   structure containing the domain partition and the quadrature rule for the boundary, 
%           since it must contain the normal vector (see msh_cartesian/msh_eval_boundary_side)
%   coeff: vector-valued function f, evaluated at the quadrature points
%
% OUTPUT:
%
%   A: assembled matrix
% 
% Copyright (C) 2015, Rafael, 2017 Vazquez
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

function varargout = op_u_otimes_n_v_otimes_n (spu, spv, msh, mshv, coeff)

  shpu = reshape (spu.shape_functions, spu.ncomp, msh.nqn, spu.nsh_max, msh.nel);
  shpv = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel);
  
  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  jacdet_weights = msh.jacdet .* msh.quad_weights .* coeff;
  
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      shpu_iel = reshape (shpu(:, :, :, iel), spu.ncomp, msh.nqn, spu.nsh_max);
      shpv_iel = reshape (shpv(:, :, :, iel), spv.ncomp, msh.nqn, spv.nsh_max);
      
      u_otimes_n_iel = zeros (spu.ncomp, spu.ncomp, msh.nqn, spu.nsh_max);
      v_otimes_n_iel = zeros (spv.ncomp, spv.ncomp, msh.nqn, spv.nsh_max);
% I need the normal from both sides
      normalu_iel = msh.normal (:, :, iel);
      normalv_iel = mshv.normal (:, :, iel);
      for ii = 1:spu.ncomp
        for jj = 1:spu.ncomp
          u_otimes_n_iel(ii,jj,:,:) = bsxfun (@times, shpu_iel(ii,:,:), normalu_iel(jj,:));
          v_otimes_n_iel(ii,jj,:,:) = bsxfun (@times, shpv_iel(ii,:,:), normalv_iel(jj,:));
        end
      end
% Should I permute it, before reshaping?
      u_otimes_n_iel = reshape (u_otimes_n_iel, spu.ncomp*spu.ncomp, msh.nqn, 1, spu.nsh_max);
      v_otimes_n_iel = reshape (v_otimes_n_iel, spv.ncomp*spv.ncomp, msh.nqn, spv.nsh_max, 1);

      jacdet_iel = reshape (jacdet_weights(:,iel), [1, msh.nqn, 1, 1]);
      v_oxn_times_jw = bsxfun (@times, jacdet_iel, v_otimes_n_iel);
      tmp1 = sum (bsxfun (@times, v_oxn_times_jw, u_otimes_n_iel), 1);
      elementary_values = reshape (sum (tmp1, 2), spv.nsh_max, spu.nsh_max);

      [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spu.connectivity(:,iel));
      indices = rows_loc & cols_loc;
      rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = rows_loc(indices);
      cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = cols_loc(indices);
      values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = elementary_values(indices);
      ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel);

    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_u_otimes_n_v_otimes_n: singular map in element number %d', iel)
    end
  end

  if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), spv.ndof, spu.ndof);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_u_otimes_n_v_otimes_n: wrong number of output arguments')
  end

end
