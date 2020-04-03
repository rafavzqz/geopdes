% OP_GRADSYMU_V_OTIMES_N: assemble the matrix A = [A(i,j)], A(i,j) = (epsilon (gradsym u)_j, (v \otimes n)_i ).
%
%   A = op_gradsymu_v_otimes_n (spu, spv, msh, coeff);
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
% Copyright (C) 2015, 2017, 2020 Rafael Vazquez
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

function varargout = op_gradsymu_v_otimes_n (spu, spv, msh, coeff)

  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], ...
		   msh.nqn, spu.nsh_max, msh.nel);

  shpv = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel);
  
  ncomp = size (gradu, 1);
  ndir = size (gradu, 2);

  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  jacdet_weights = msh.jacdet .* msh.quad_weights .* coeff;
  
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      gradu_iel = reshape (gradu(:,:,:,:,iel), spu.ncomp, ndir, msh.nqn, spu.nsh_max);
      gradu_iel = 0.5 * (gradu_iel + permute(gradu_iel, [2 1 3 4]));
      gradu_iel = reshape (gradu_iel, spu.ncomp*ndir, msh.nqn, 1, spu.nsh_max);
      shpv_iel = reshape (shpv(:, :, :, iel), spv.ncomp, msh.nqn, spv.nsh_max);
      
      v_otimes_n_iel = zeros (ncomp, ndir, msh.nqn, spv.nsh_max);
      normal_iel = msh.normal (:, :, iel);
      for ii = 1:spv.ncomp
        for jj = 1:ndir
          v_otimes_n_iel(ii,jj,:,:) = bsxfun (@times, shpv_iel(ii,:,:), normal_iel(jj,:));
        end
      end
      v_otimes_n_iel = reshape (v_otimes_n_iel, ncomp*ndir, msh.nqn, spv.nsh_max, 1);
        
      jacdet_iel = reshape (jacdet_weights(:,iel), [1, msh.nqn, 1]);
      v_oxn_times_jw = bsxfun (@times, jacdet_iel, v_otimes_n_iel);
      tmp1 = sum (bsxfun (@times, v_oxn_times_jw, gradu_iel), 1);
      elementary_values = reshape (sum (tmp1, 2), spv.nsh_max, spu.nsh_max);

      [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spu.connectivity(:,iel));
      indices = rows_loc & cols_loc;
      rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = rows_loc(indices);
      cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = cols_loc(indices);
      values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = elementary_values(indices);
      ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel);
      
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradsymu_v_otimes_n: singular map in element number %d', iel)
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
    error ('op_gradsymu_v_otimes_n: wrong number of output arguments')
  end

end
