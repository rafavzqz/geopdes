% OP_SU_EV: assemble the matrix A = [a(i,j)], a(i,j) = 1/2 (sigma (u_j), epsilon (v_i)).
%
%   mat = op_su_ev (spu, spv, msh, lambda, mu);
%   [rows, cols, values] = op_su_ev (spu, spv, msh, lambda, mu);
%
% INPUT:
%    
%   spu:   structure representing the space of trial functions (see sp_bspline_2d/sp_evaluate_col)
%   spv:   structure representing the space of test functions (see sp_bspline_2d/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_2d/msh_evaluate_col)
%   lambda, mu: Lame' coefficients evaluated at the quadrature points
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011 Rafael Vazquez
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

function varargout = op_su_ev (spu, spv, msh, lambda, mu)

  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);

  ndir = size (gradu, 2);

  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      jacdet_weights = msh.jacdet(:, iel) .* msh.quad_weights(:, iel);
      jacdet_weights_mu = reshape (jacdet_weights .* mu(:, iel), 1, msh.nqn);
      jacdet_weights_lambda = reshape (jacdet_weights .* lambda(:, iel), msh.nqn, 1);
      
      gradu_iel = permute (gradu(:, :, :, 1:spu.nsh(iel), iel), [1 2 4 3]);
      epsu_iel = (gradu_iel + permute (gradu_iel, [2 1 3 4]))/2;
      epsu_iel = reshape (epsu_iel, spu.ncomp * ndir, spu.nsh(iel), msh.nqn);
      epsu_iel = permute (epsu_iel, [1 3 2]);

      gradv_iel = permute (gradv(:, :, :, 1:spv.nsh(iel), iel), [1 2 4 3]);
      epsv_iel = (gradv_iel + permute (gradv_iel, [2 1 3 4]))/2;
      epsv_iel = reshape (epsv_iel, spv.ncomp * ndir, spv.nsh(iel), msh.nqn);
      epsv_iel = permute (epsv_iel, [1 3 2]);

      divv_iel = spv.shape_function_divs(:,1:spv.nsh(iel),iel);
      divu_iel = spu.shape_function_divs(:,1:spu.nsh(iel),iel);

      epsv_times_jw = bsxfun (@times, jacdet_weights_mu, epsv_iel);
      divv_times_jw = bsxfun (@times, jacdet_weights_lambda, divv_iel);
      for idof = 1:spv.nsh(iel)
        rows(ncounter+(1:spu.nsh(iel))) = spv.connectivity(idof, iel);
        cols(ncounter+(1:spu.nsh(iel))) = spu.connectivity(1:spu.nsh(iel), iel);

        aux_val1 = bsxfun (@times, epsv_times_jw(:,:,idof), epsu_iel);
        aux_val2 = bsxfun (@times, divv_times_jw(:,idof), divu_iel);
        values(ncounter+(1:spu.nsh(iel))) = 2 * sum (sum (aux_val1, 2), 1);
        values(ncounter+(1:spu.nsh(iel))) = values(ncounter+(1:spu.nsh(iel))) + sum (aux_val2, 1).';

        ncounter = ncounter + spu.nsh(iel);
      end
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_su_ev: singular map in element number %d', iel)
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
    error ('op_su_ev: wrong number of output arguments')
  end

end

%% COPY OF THE FIRST VERSION OF THE FUNCTION (MORE UNDERSTANDABLE)
% 
% function mat = op_su_ev (spu, spv, msh, lambda, mu)
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
%         ishg  = gradv(:,:,:,idof,iel);
%         ishgt = permute (ishg, [2, 1, 3]);
%         ieps  = reshape(ishg + ishgt, spv.ncomp * ndir, [])/2;
%         idiv  = spv.shape_function_divs(:, idof, iel);
%         for jdof = 1:spu.nsh(iel) 
%           jshg  = gradu(:,:,:,jdof,iel);
%           jshgt = permute (jshg, [2, 1, 3]);
%           jeps  = reshape(jshg + jshgt, spu.ncomp * ndir, [])/2;
%           jdiv  = spu.shape_function_divs(:, jdof, iel);
%  % The cycle on the quadrature points is vectorized         
%           mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
%               sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* ...
%                    (2 * sum (ieps .* jeps, 1).' .* mu(:,iel)  + ...
%                     (idiv .* jdiv) .* lambda(:,iel)));
%         end
%       end
%       mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
%         mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) + mat_loc;
%     else
%       warning ('geopdes:jacdet_zero_at_quad_node', 'op_su_ev: singular map in element number %d', iel)
%     end
%   end
% 
% end
