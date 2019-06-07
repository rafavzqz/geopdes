% OP_UGRADU_JAC_V: assemble the second term of the jacobian of the
% convective term (see solve_navier_stokes for usage)
%
%   mat = op_ugradu_jac_v (spu, spv, msh, coeff);
%   [rows, cols, values] = op_ugradu_jac (spu, spv, msh, uder);
%
% INPUT:
%
%  spu:   structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%  spv:   structure representing the space of test functions  (see sp_scalar/sp_evaluate_col)
%  msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%  uder:  derivative of the velocity, evaluated in all the degrees of
%  freedom
%
% OUTPUT:
%
%  mat:    assembled jacobian matrix
%  rows:   row indices of the nonzero entries
%  cols:   column indices of the nonzero entries
%  values: values of the nonzero entries
% 
% Copyright (C) 2018 Luca Coradello, Luca Pegolotti
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

function varargout = op_ugradu_jac_v (spu, spv, msh, uder)

  gradu_shpu = sum(bsxfun( @times, reshape(spu.shape_functions, ...
                   [1 size(spu.shape_functions)]),uder),2);

  shpu = reshape (gradu_shpu, spu.ncomp, msh.nqn, spu.nsh_max, msh.nel);
  shpv = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel);
  
  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  jacdet_weights = msh.jacdet .* msh.quad_weights;
  
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
      shpu_iel = reshape (shpu(:, :, 1:spu.nsh(iel), iel), spu.ncomp, msh.nqn, 1, spu.nsh(iel));
      shpv_iel = reshape (shpv(:, :, 1:spv.nsh(iel), iel), spv.ncomp, msh.nqn, spv.nsh(iel), 1);

      jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);

      jacdet_shpu = bsxfun (@times, jacdet_iel, shpu_iel);
      tmp1 = sum (bsxfun (@times, jacdet_shpu, shpv_iel), 1);
      values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = reshape (sum (tmp1, 2), spv.nsh(iel), spu.nsh(iel));

      [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spu.connectivity(:,iel));
      rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = rows_loc;
      cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = cols_loc;
      ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel);
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_u_v: singular map in element number %d', iel)
    end
  end

  if (nargout == 1 || nargout == 0)
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