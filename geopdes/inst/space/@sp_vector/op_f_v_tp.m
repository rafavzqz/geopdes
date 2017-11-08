% OP_F_V_TP: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i), exploiting the tensor product structure.
%
%   rhs = op_f_v_tp (spv, msh, coeff);
%
% INPUT:
%     
%   spv:   object representing the function space (see sp_vector)
%   msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%   coeff: function handle to compute the source function
%
% OUTPUT:
%
%   rhs: assembled right-hand side
% 
% Copyright (C) 2011, 2017 Rafael Vazquez
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

function rhs = op_f_v_tp (space, msh, coeff)

  for icomp = 1:space.ncomp_param
    for idim = 1:msh.ndim
      size1 = size (space.scalar_spaces{icomp}.sp_univ(idim).connectivity);
      if (size1(2) ~= msh.nel_dir(idim))
        error ('The discrete space is not associated to the mesh')
      end
    end
  end

  rhs = zeros (space.ndof, 1);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col);

    for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end

    rhs = rhs + op_f_v (sp_col, msh_col, coeff (x{:}));
  end

end
