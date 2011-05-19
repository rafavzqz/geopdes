% OP_U_V: assemble the mass matrix M = [m(i,j)], m(i,j) = (mu u_j, v_i).
%
%   mat = op_u_v (spu, spv, msh, coeff);
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
%  mat: assembled mass matrix
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

function mat = op_u_v (spu, spv, msh, coeff)
  
  shpu = reshape (spu.shape_functions, spu.ncomp, msh.nqn, spu.nsh_max, msh.nel);
  shpv = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel);
  
  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
      jacdet_weights = msh.jacdet(:, iel) .* ...
                       msh.quad_weights(:, iel) .* coeff(:, iel);

      shpv_iel = reshape (shpv(:, :, :, iel), spv.ncomp, msh.nqn, spv.nsh_max);
      shpu_iel = reshape (shpu(:, :, :, iel), spu.ncomp, msh.nqn, spu.nsh_max);
      for idof = 1:spv.nsh(iel)
        ishp = shpv_iel(:, :, idof);
        for jdof = 1:spu.nsh(iel)
          ncounter = ncounter + 1;
          rows(ncounter) = spv.connectivity(idof, iel);
          cols(ncounter) = spu.connectivity(jdof, iel);

          jshp = shpu_iel(:, :, jdof);

          values(ncounter) = sum (jacdet_weights .* sum (ishp .* jshp, 1).');
        end
      end
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_u_v: singular map in element number %d', iel)
    end
  end

  mat = sparse (rows, cols, values, spv.ndof, spu.ndof);

end

