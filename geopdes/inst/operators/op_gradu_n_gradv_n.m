% OP_GRADU_N_GRADV_N: assemble the matrix A = [a(i,j)], a(i,j) = (epsilon (grad u n)_j, (grad v n)_i), with n the normal vector.
%
%   mat = op_gradu_n_gradv_n (spu, spv, msh, epsilon);
%   [rows, cols, values] = op_gradu_n_gradv_n (spu, spv, msh, epsilon);
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%   spv:   structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule for the boundary, 
%           since it must contain the normal vector (see msh_cartesian/msh_eval_boundary_side)
%   epsilon: coefficient
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2014 Adriano Cortes
% Copyright (C) 2014, 2017, 2018 Rafael Vazquez
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

function varargout = op_gradu_n_gradv_n (spu, spv, msh, coeff)

  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], ...
		   msh.nqn, spu.nsh_max, msh.nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], ...
		   msh.nqn, spv.nsh_max, msh.nel);

  ndim = size (gradu, 2);

  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  jacdet_weights = msh.jacdet .* msh.quad_weights .* coeff;
  
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
      gradu_iel = gradu(:,:,:,:,iel);
      gradv_iel = gradv(:,:,:,:,iel);
      normal_iel = reshape (msh.normal(:,:,iel), [1, ndim, msh.nqn]);

      gradu_n = reshape (sum (bsxfun (@times, gradu_iel, normal_iel), 2), spu.ncomp, msh.nqn, 1, spu.nsh_max);
      gradv_n = reshape (sum (bsxfun (@times, gradv_iel, normal_iel), 2), spv.ncomp, msh.nqn, spv.nsh_max, 1);

      jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);

      gradv_n_times_jw = bsxfun (@times, jacdet_iel, gradv_n);
      tmp1 = sum (bsxfun (@times, gradv_n_times_jw, gradu_n), 1);
      elementary_values = reshape (sum (tmp1, 2), spv.nsh_max, spu.nsh_max);
      
      [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spu.connectivity(:,iel));
      indices = rows_loc & cols_loc;
      rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = rows_loc(indices);
      cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = cols_loc(indices);
      values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = elementary_values(indices);
      ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel);
      
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradv_n_u: singular map in element number %d', iel)
    end
  end
  

  if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows, cols, values, spv.ndof, spu.ndof);
  elseif (nargout == 3)
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = values;
  else
    error ('op_gradu_n_gradv_n: wrong number of output arguments')
  end
  
  
end
