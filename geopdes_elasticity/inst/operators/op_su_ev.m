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

function varargout = op_su_ev (spu, spv, msh, lambda, mu, element_list)

  if (nargin == 5)
    element_list = 1:msh.nel;
  end
  nel = numel (element_list);

  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, nel);

  ndir = size (gradu, 2);

  rows = zeros (nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (nel * spu.nsh_max * spv.nsh_max, 1);

  ncounter = 0;
  for iel = 1:nel
    iel_glob = element_list (iel);
    if (all (msh.jacdet(:, iel_glob)))
      jacdet_weights = msh.jacdet(:, iel_glob) .* msh.quad_weights(:, iel_glob);
      jacdet_weights_mu = repmat (jacdet_weights .* mu(:, iel_glob), [1,spu.nsh(iel)]);
      jacdet_weights_lambda = repmat (jacdet_weights .* lambda(:, iel_glob), [1,spu.nsh(iel)]);

      gradu_iel = permute (gradu(:, :, :, 1:spu.nsh(iel), iel), [1 2 4 3]);
      epsu_iel = (gradu_iel + permute (gradu_iel, [2 1 3 4]))/2;
      epsu_iel = reshape (epsu_iel, spu.ncomp * ndir, spu.nsh(iel), msh.nqn);
      epsu_iel = permute (epsu_iel, [1 3 2]);

      gradv_iel = permute (gradv(:, :, :, 1:spv.nsh(iel), iel), [1 2 4 3]);
      epsv_iel = (gradv_iel + permute (gradv_iel, [2 1 3 4]))/2;
      epsv_iel = reshape (epsv_iel, spv.ncomp * ndir, spv.nsh(iel), msh.nqn);
      epsv_iel = permute (epsv_iel, [1 3 2]);

      divu_iel = spu.shape_function_divs(:,:,iel);

      for idof = 1:spv.nsh(iel)
        ieps  = repmat (epsv_iel(:, :, idof), [1 1 spu.nsh(iel)]);;
        idiv  = repmat (spv.shape_function_divs(:, idof, iel), [1 spu.nsh(iel)]);

        rows(ncounter+(1:spu.nsh(iel))) = spv.connectivity(idof, iel);
        cols(ncounter+(1:spu.nsh(iel))) = spu.connectivity(1:spu.nsh(iel), iel);

        val1 = 2 * sum (jacdet_weights_mu .* ...
          reshape (sum (ieps .* epsu_iel, 1), msh.nqn, spu.nsh(iel)), 1);
        val2 = sum (jacdet_weights_lambda .* idiv .* divu_iel, 1);

        values(ncounter+(1:spu.nsh(iel))) = val1 + val2;
        ncounter = ncounter + spu.nsh(iel);
      end
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_su_ev: singular map in element number %d', iel)
    end
  end

  if (nargout == 1)
    varargout{1} = sparse (rows, cols, values, spv.ndof, spu.ndof);
  elseif (nargout == 3)
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = values;
  else
    error ('op_su_ev: wrong number of output arguments')
  end

end
