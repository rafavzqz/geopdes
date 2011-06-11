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
%   mat:    assembled stiffness matrix
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

function varargout = op_gradu_gradv (spu, spv, msh, coeff)
  
  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);

  ndir = size (gradu, 2);

  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      jacdet_weights = msh.jacdet(:, iel) .* ...
                       msh.quad_weights(:, iel) .* coeff(:, iel);

      gradu_iel = permute (gradu(:, :, :, :, iel), [1 2 4 3]);
      gradu_iel = reshape (gradu_iel, spu.ncomp * ndir, spu.nsh_max, []);
      gradu_iel = permute (gradu_iel, [1 3 2]);

      gradv_iel = permute (gradv(:, :, :, :, iel), [1 2 4 3]);
      gradv_iel = reshape (gradv_iel, spv.ncomp * ndir, spv.nsh_max, []);
      gradv_iel = permute (gradv_iel, [1 3 2]);

      for idof = 1:spv.nsh(iel)
        ishg = gradv_iel(:, :, idof);
        for jdof = 1:spu.nsh(iel) 
          ncounter = ncounter + 1;
          rows(ncounter) = spv.connectivity(idof, iel);
          cols(ncounter) = spu.connectivity(jdof, iel);

          jshg = gradu_iel(:, :, jdof);

          values(ncounter) = sum (jacdet_weights .* sum (ishg .* jshg, 1).');
        end
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