% OP_VEL_DOT_GRADU_V: assemble the matrix B = [b(i,j)], b(i,j) = (V * grad u_j, v_i).
%  The current version only works for the scalar case.
%
%   mat = op_vel_dot_gradu_v (spu, spv, msh, Vel);
%   [rows, cols, values] = op_vel_dot_gradu_v (spu, spv, msh, Vel);
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%   spv:   structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   Vel:   advection field
%
% OUTPUT:
%
%   mat:    assemble advection matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2014, 2017 Rafael Vazquez
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

function varargout = op_vel_dot_gradu_v (spu, spv, msh, coeff)

  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], ...
		   msh.nqn, spu.nsh_max, msh.nel);
  shpv = reshape (spv.shape_functions, spv.ncomp, ...
                  msh.nqn, spv.nsh_max, msh.nel);
  
  ndir = size (gradu, 2);

  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  jacdet_weights = msh.jacdet .* msh.quad_weights;
  
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      gradu_iel = reshape (gradu(:,:,:,:,iel), spu.ncomp*ndir, msh.nqn, 1, spu.nsh_max);

      coeff_iel = reshape (coeff(:, :, iel), [ndir, msh.nqn, 1, 1]);

      vel_dot_gradu_iel = bsxfun (@times, gradu_iel, coeff_iel);

      shpv_iel = reshape (shpv(:, :, :, iel), spv.ncomp, msh.nqn, spv.nsh_max, 1);
      jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);

      jacdet_shpv = bsxfun (@times, jacdet_iel, shpv_iel);
      tmp1 = bsxfun (@times, jacdet_shpv, vel_dot_gradu_iel);
      elementary_values = reshape (sum (sum (tmp1, 1), 2), spv.nsh_max, spu.nsh_max);

      [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spu.connectivity(:,iel));
      indices = rows_loc & cols_loc;
      rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = rows_loc(indices);
      cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = cols_loc(indices);
      values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = elementary_values(indices);
      ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel);

    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_vel_dot_gradu_v: singular map in element number %d', iel)
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
    error ('op_vel_dot_gradu_v: wrong number of output arguments')
  end

end
