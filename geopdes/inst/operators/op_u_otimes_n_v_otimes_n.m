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
% Copyright (C) 2015, Rafael Vazquez
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

% I don't need the theta-theta component of the gradient
  shpu = reshape (spu.shape_functions, spu.ncomp, msh.nqn, spu.nsh_max, msh.nel);
  shpv = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel);
  
  ncomp = spu.ncomp;

  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      jacdet_weights = reshape (msh.jacdet(:, iel) .* ...
                       msh.quad_weights(:, iel) .* coeff(:, iel), 1, msh.nqn);

% I cannot use spu.ncomp here

      shpu_iel = reshape (shpu(:, :, 1:spu.nsh(iel), iel), spu.ncomp, msh.nqn, spu.nsh(iel));
      shpv_iel = reshape (shpv(:, :, 1:spv.nsh(iel), iel), spv.ncomp, msh.nqn, spv.nsh(iel));
      u_otimes_n_iel = zeros (ncomp, ncomp, msh.nqn, spu.nsh(iel));
      v_otimes_n_iel = zeros (ncomp, ncomp, msh.nqn, spv.nsh(iel));
% I need the normal from both sides
      normalu = msh.normal (:, :, iel);
      normalv = mshv.normal (:, :, iel);
      for ii = 1:ncomp
        for jj = 1:ncomp
          u_otimes_n_iel(ii,jj,:,:) = bsxfun (@times, shpu_iel(ii,:,:), normalu(jj,:));
          v_otimes_n_iel(ii,jj,:,:) = bsxfun (@times, shpv_iel(ii,:,:), normalv(jj,:));
        end
      end
% Should I permute it, before reshaping?
      u_otimes_n_iel = reshape (u_otimes_n_iel, ncomp*ncomp, msh.nqn, spv.nsh(iel));
      v_otimes_n_iel = reshape (v_otimes_n_iel, ncomp*ncomp, msh.nqn, spv.nsh(iel));

      v_oxn_times_jw = bsxfun (@times, jacdet_weights, v_otimes_n_iel);
      for idof = 1:spv.nsh(iel)
        rows(ncounter+(1:spu.nsh(iel))) = spv.connectivity(idof, iel);
        cols(ncounter+(1:spu.nsh(iel))) = spu.connectivity(1:spu.nsh(iel), iel);

        aux_val = bsxfun (@times, v_oxn_times_jw(:,:,idof), u_otimes_n_iel);
        values(ncounter+(1:spu.nsh(iel))) = sum (sum (aux_val, 2), 1);
        ncounter = ncounter + spu.nsh(iel);
      end
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradu_gradv: singular map in element number %d', iel)
    end
  end

  if (nargout == 1)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), spv.ndof, spu.ndof);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_gradu_gradv: wrong number of output arguments')
  end

end



%% COPY OF THE FIRST VERSION OF THE FUNCTION (MORE UNDERSTANDABLE)
% 
% function mat = op_gradu_gradv (spu, spv, msh, coeff)
%   
%   mat = spalloc (spv.ndof, spu.ndof, 1);
%   
%   gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
%   gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);
% 
%   ndir = size (gradu, 2);
% 
%   for iel = 1:msh.nel
%     if (all (msh.jacdet(:,iel)))
%       mat_loc = zeros (spv.nsh(iel), spu.nsh(iel));
%       for idof = 1:spv.nsh(iel)
%         ishg = reshape(gradv(:,:,:,idof,iel),spv.ncomp * ndir, []);
%         for jdof = 1:spu.nsh(iel) 
%           jshg = reshape(gradu(:,:,:,jdof,iel),spu.ncomp * ndir, []);
% % The cycle on the quadrature points is vectorized
%           %for inode = 1:msh.nqn
%           mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
%              sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* ...
%                   sum (ishg .* jshg, 1).' .* coeff(:,iel));
%           %end  
%         end
%       end
%       mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
%         mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) + mat_loc;
%     else
%       warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradu_gradv: singular map in element number %d', iel)
%     end
%   end
% 
% end
