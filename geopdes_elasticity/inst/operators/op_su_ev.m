% OP_SU_EV: assemble the matrix A = [a(i,j)], a(i,j) = 1/2 (sigma (u_j), epsilon (v_i)).
%
%   mat = op_su_ev (spu, spv, msh, lambda, mu);
%
% INPUT:
%    
%   spu:        structure representing the space of trial functions (see sp_bspline_2d_phys)
%   spv:        structure representing the space of test functions  (see sp_bspline_2d_phys)
%   msh:        structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%   lambda, mu: Lame' coefficients
%
% OUTPUT:
%
%   mat: assembled matrix
% 
% Copyright (C) 2009, 2010 Carlo de Falco
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

function mat = op_su_ev (spu, spv, msh, lambda, mu)
  
  mat = spalloc (spv.ndof, spu.ndof, 1);
  
  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);

  ndir = size (gradu, 2);

  for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
      mat_loc = zeros (spv.nsh(iel), spu.nsh(iel));
      for idof = 1:spv.nsh(iel)
        ishg  = gradv(:,:,:,idof,iel);
        ishgt = permute (ishg, [2, 1, 3]);
        ieps  = reshape(ishg + ishgt, spv.ncomp * ndir, [])/2;
        idiv  = spv.shape_function_divs(:, idof, iel);
        for jdof = 1:spu.nsh(iel) 
          jshg  = gradu(:,:,:,jdof,iel);
          jshgt = permute (jshg, [2, 1, 3]);
          jeps  = reshape(jshg + jshgt, spu.ncomp * ndir, [])/2;
          jdiv  = spu.shape_function_divs(:, jdof, iel);
          
          mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
              sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* ...
                   (2 * sum (ieps .* jeps, 1).' .* mu(:,iel)  + ...
                    (idiv .* jdiv) .* lambda(:,iel)));
        end
      end
      mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
        mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) + mat_loc;
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_su_ev: singular map in element number %d', iel)
    end
  end

end
