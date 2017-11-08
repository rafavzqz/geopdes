% OP_UXN_VXN_3D: assemble the matrix M = [m(i,j)], m(i,j) = (mu u_j x n, v_i x n), with n the exterior normal vector.
%
%   mat = op_uxn_vxn_3d (spu, spv, msh, coeff);
%   [rows, cols, values] = op_uxn_vxn_3d (spu, spv, msh, coeff);
%
% INPUT:
%   
%  spu:   structure representing the space of trial functions (see sp_vector/sp_eval_boundary_side)
%  spv:   structure representing the space of test functions  (see sp_vector/sp_eval_boundary_side)
%  msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_eval_boundary_side)
%  coeff: physical parameter
%
% OUTPUT:
%
%  mat:    assembled matrix
%  rows:   row indices of the nonzero entries
%  cols:   column indices of the nonzero entries
%  values: values of the nonzero entries
% 
% Copyright (C) 2009, 2010 Carlo de Falco, Rafael Vazquez
% Copyright (C) 2011, 2017 Rafael Vazquez
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

function varargout = op_uxn_vxn_3d (spu, spv, msh, coeff)
  
  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  jacdet_weights = msh.jacdet .* msh.quad_weights .* coeff;

  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      shpu_iel = reshape (spu.shape_functions(:, :, :, iel), spu.ncomp, msh.nqn, 1, spu.nsh_max);
      shpv_iel = reshape (spv.shape_functions(:, :, :, iel), spv.ncomp, msh.nqn, spv.nsh_max, 1);

      normal_iel = reshape (msh.normal(:,:,iel), 3, msh.nqn);

      shpu_x_n = bsxfun (@cross, shpu_iel, normal_iel);
      shpv_x_n = bsxfun (@cross, shpv_iel, normal_iel);
      
      jacdet_iel = reshape (jacdet_weights(:,iel), [1, msh.nqn, 1]);
      jacdet_shpu = bsxfun (@times, jacdet_iel, shpu_x_n);
      tmp1 = sum (bsxfun (@times, jacdet_shpu, shpv_x_n), 1);
      elementary_values = reshape (sum (tmp1, 2), spv.nsh_max, spu.nsh_max);
      
      [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spu.connectivity(:,iel));
      indices = rows_loc & cols_loc;
      rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = rows_loc(indices);
      cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = cols_loc(indices);
      values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = elementary_values(indices);
      ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel);

    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_uxn_vxn_3d: singular map in element number %d', iel)
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
    error ('op_uxn_vxn_3d: wrong number of output arguments')
  end

end

