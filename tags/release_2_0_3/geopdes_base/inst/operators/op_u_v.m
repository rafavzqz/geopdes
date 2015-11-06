% OP_U_V: assemble the mass matrix M = [m(i,j)], m(i,j) = (mu u_j, v_i).
%
%   mat = op_u_v (spu, spv, msh, coeff);
%   [rows, cols, values] = op_u_v (spu, spv, msh, coeff);
%
% INPUT:
%
%  spu:   structure representing the space of trial functions (see sp_bspline_2d_phys)
%  spv:   structure representing the space of test functions  (see sp_bspline_2d_phys)
%  msh:   structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%  coeff: reaction coefficient
%
% OUTPUT:
%
%  mat:    assembled mass matrix
%  rows:   row indices of the nonzero entries
%  cols:   column indices of the nonzero entries
%  values: values of the nonzero entries
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

function varargout = op_u_v (spu, spv, msh, coeff)
  
  shpu = reshape (spu.shape_functions, spu.ncomp, msh.nqn, spu.nsh_max, msh.nel);
  shpv = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel);
  
  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
      jacdet_weights = reshape (msh.jacdet(:, iel) .* ...
                       msh.quad_weights(:, iel) .* coeff(:, iel), 1, msh.nqn);

      shpv_iel = reshape (shpv(:, :, 1:spv.nsh(iel), iel), spv.ncomp, msh.nqn, spv.nsh(iel));
      shpu_iel = reshape (shpu(:, :, 1:spu.nsh(iel), iel), spu.ncomp, msh.nqn, spu.nsh(iel));

      shpv_times_jw = bsxfun (@times, jacdet_weights, shpv_iel);
      for idof = 1:spv.nsh(iel)
        rows(ncounter+(1:spu.nsh(iel))) = spv.connectivity(idof, iel);
        cols(ncounter+(1:spu.nsh(iel))) = spu.connectivity(1:spu.nsh(iel), iel);

        aux_val = bsxfun (@times, shpv_times_jw(:,:,idof), shpu_iel);
        values(ncounter+(1:spu.nsh(iel))) = sum (sum (aux_val, 2), 1);
        ncounter = ncounter + spu.nsh(iel);
      end
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_u_v: singular map in element number %d', iel)
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
    error ('op_u_v: wrong number of output arguments')
  end

end


%% COPY OF THE FIRST VERSION OF THE FUNCTION (MORE UNDERSTANDABLE)
% 
% function mat = op_u_v (spu, spv, msh, coeff)
%   
%   mat = spalloc(spv.ndof, spu.ndof, 1);
%   shpu = reshape (spu.shape_functions, spu.ncomp, msh.nqn, spu.nsh_max, msh.nel);
%   shpv = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel);
%   
%   for iel = 1:msh.nel
%     if (all (msh.jacdet(:,iel)))
%       mat_loc = zeros (spv.nsh(iel), spu.nsh(iel));
%       for idof = 1:spv.nsh(iel)
%         ishp = reshape (shpv(:,:,idof,iel), spv.ncomp, []);
%         for jdof = 1:spu.nsh(iel)
%           jshp = reshape (shpu(:,:,jdof,iel), spu.ncomp, []);
% % The cycle on the quadrature points is vectorized
%           %for inode = 1:msh.nqn
%             mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
%               sum (msh.jacdet(:, iel) .* msh.quad_weights(:, iel) .* ...
%               sum(ishp .* jshp, 1).' .* coeff(:, iel));
%           %end  
%         end
%       end
%       mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
%         mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) +  mat_loc;
%     else
%       warning ('geopdes:jacdet_zero_at_quad_node', 'op_u_v: singular map in element number %d', iel)
%     end
%   end
% 
% end
