% OP_GRADVN_U: assemble the matrix A = [a(i,j)], a(i,j) = (epsilon (grad v n)_j, u_i), with n the normal vector.
%
%   mat = op_gradvn_u (spu, spv, msh, epsilon);
%   [rows, cols, values] = op_gradu_gradv (spu, spv, msh, epsilon);
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_bspline_2d/sp_evaluate_col)
%   spv:   structure representing the space of test functions (see sp_bspline_2d/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule for the boundary, 
%           since it must contain the normal vector (see msh_2d/msh_eval_boundary_side)
%   epsilon: coefficient
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2014 Adriano Cortes, Rafael Vazquez
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

function varargout = op_gradv_n_u (spu, spv, msh, coeff)

  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], ...
		   msh.nqn, spv.nsh_max, msh.nel);

  ndim = size (gradv, 2);

  shpu = reshape (spu.shape_functions, spu.ncomp, msh.nqn, spu.nsh_max, msh.nel);

  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))

      jacdet_weights = reshape (msh.jacdet(:, iel) .* ...
                       msh.quad_weights(:, iel) .* coeff(:, iel), 1, msh.nqn);

      gradv_iel = gradv(:, :, :, 1:spv.nsh(iel), iel);

      normal_iel = reshape(msh.normal(:,:,iel), [1 ndim msh.nqn]);

      shpu_iel = reshape (shpu(:, :, 1:spu.nsh(iel), iel), spu.ncomp, msh.nqn, spu.nsh(iel));

      gradv_n = reshape (sum (bsxfun(@times, gradv_iel, normal_iel), 2), ...
                        [spv.ncomp, msh.nqn, spv.nsh(iel)]);

      gradv_n_times_jw = bsxfun (@times, jacdet_weights, gradv_n);
      
      for idof = 1:spv.nsh(iel)
        rows(ncounter+(1:spu.nsh(iel))) = spv.connectivity(idof, iel);
        cols(ncounter+(1:spu.nsh(iel))) = spu.connectivity(1:spu.nsh(iel), iel);

        aux_val = bsxfun (@times, gradv_n_times_jw(:,:,idof), shpu_iel);
        values(ncounter+(1:spu.nsh(iel))) = sum (sum (aux_val, 2), 1);
        ncounter = ncounter + spu.nsh(iel);
      end

    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradv_n_u: singular map in element number %d', iel)
    end
  end
  

  if (nargout == 1)
    varargout{1} = sparse (rows, cols, values, spv.ndof, spu.ndof);
  elseif (nargout == 3)
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = values;
  else
    error ('op_gradv_n_u: wrong number of output arguments')
  end
  
  
end
