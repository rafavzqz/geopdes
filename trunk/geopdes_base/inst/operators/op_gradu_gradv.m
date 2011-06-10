% OP_GRADU_GRADV: assemble the stiffness matrix A = [a(i,j)], a(i,j) = (epsilon grad u_j, grad v_i).
%
%   mat = op_gradu_gradv (spu, spv, msh, epsilon);
%   [rows, cols, values] = op_gradu_gradv (spu, spv, msh, epsilon);
%
% INPUT:
%    
%   spu:     structure representing the space of trial functions (see sp_bspline_2d_phys)
%   spv:     structure representing the space of test functions  (see sp_bspline_2d_phys)
%   msh:     structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%   epsilon: diffusion coefficient
%
% OUTPUT:
%
%   mat: assembled stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, Rafael Vazquez
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

function varargout = op_gradu_gradv (spu, spv, msh, coeff, element_list)

  if (nargin == 4)
    element_list = 1:msh.nel;
  end
  nel = numel (element_list);

  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], ...
		   msh.nqn, spu.nsh_max, nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], ...
		   msh.nqn, spv.nsh_max, nel);

  ndir = size (gradu, 2);

  rows = zeros (nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (nel * spu.nsh_max * spv.nsh_max, 1);

  ncounter = 0;
  for iel = 1:nel
    iel_glob = element_list (iel);
    if (all (msh.jacdet(:, iel_glob)))
      jacdet_weights = msh.jacdet(:, iel_glob) .* ...
                       msh.quad_weights(:, iel_glob) .* coeff(:, iel_glob);
      jacdet_weights = repmat (jacdet_weights, [1,spu.nsh(iel)]);

      gradu_iel = permute (gradu(:, :, :, 1:spu.nsh(iel), iel), [1 2 4 3]);
      gradu_iel = reshape (gradu_iel, spu.ncomp * ndir, spu.nsh(iel), msh.nqn);
      gradu_iel = permute (gradu_iel, [1 3 2]);

      gradv_iel = permute (gradv(:, :, :, 1:spv.nsh(iel), iel), [1 2 4 3]);
      gradv_iel = reshape (gradv_iel, spv.ncomp * ndir, spv.nsh(iel), msh.nqn);
      gradv_iel = permute (gradv_iel, [1 3 2]);

      for idof = 1:spv.nsh(iel)
        ishg = repmat (gradv_iel(:, :, idof), [1 1 spu.nsh(iel)]);

        rows(ncounter+(1:spu.nsh(iel))) = spv.connectivity(idof, iel);
        cols(ncounter+(1:spu.nsh(iel))) = spu.connectivity(1:spu.nsh(iel), iel);

        values(ncounter+(1:spu.nsh(iel))) = sum (jacdet_weights .* ...
          reshape (sum (ishg .* gradu_iel, 1), msh.nqn, spu.nsh(iel)), 1);
        ncounter = ncounter + spu.nsh(iel);
      end
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradu_gradv: singular map in element number %d', iel)
    end
  end

  if (nargout == 1)
    varargout{1} = sparse (rows, cols, values, spv.ndof, spu.ndof);
  elseif (nargout == 3)
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = values;
  else
    error ('op_gradu_gradv: wrong number of output arguments')
  end

end
