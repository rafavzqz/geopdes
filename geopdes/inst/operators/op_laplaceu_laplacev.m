% OP_LAPLACEU_LAPLACEV: assemble the matrix A = [a(i,j)], a(i,j) = (epsilon laplace u_j, laplace v_i).
%
%   mat = op_laplaceu_laplacev (spu, spv, msh, epsilon);
%   [rows, cols, values] = op_laplaceu_laplacev (spu, spv, msh, epsilon);
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%   spv:   structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   epsilon: diffusion coefficient
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, Rafael Vazquez
% Copyright (C) 2013, Marco Pingaro
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

function mat = op_laplaceu_laplacev (spu, spv, msh, coeff)
   mat = spalloc (spv.ndof, spu.ndof, 1);
   
   laplaceu = reshape (spu.shape_function_hessians, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
   laplacev = reshape (spv.shape_function_hessians, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);
   
   ndir = size (laplaceu, 2);
 
   for iel = 1:msh.nel
     if (all (msh.jacdet(:,iel)))
       mat_loc = zeros (spv.nsh(iel), spu.nsh(iel));
       for idof = 1:spv.nsh(iel)
         ishg = reshape(laplacev(:,:,:,idof,iel),spv.ncomp * ndir, []);
         for jdof = 1:spu.nsh(iel) 
           jshg = reshape(laplaceu(:,:,:,jdof,iel),spu.ncomp * ndir, []);
           % The cycle on the quadrature points is vectorized
           %for inode = 1:msh.nqn
           mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
             sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* ...
             sum (ishg([1,1,4,4],:) .* jshg([1,4,1,4],:), 1).' .* coeff(:,iel));
           %end
         end
       end
       mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
         mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) + mat_loc;
     else
       warning ('geopdes:jacdet_zero_at_quad_node',...
           'op_laplaceu_laplacev: singular map in element number %d', iel)
     end
   end
end
