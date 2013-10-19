% OP_F_V_TP: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i), exploiting the tensor product structure.
%
%   rhs = op_f_v_tp (spv, msh, coeff);
%   rhs = op_f_v_tp (spv, msh, coeff_sp, coeff_dofs);
%
% INPUT:
%     
%   spv:        class representing the function space (see sp_bspline_2d)
%   msh:        class defining the domain partition and the quadrature rule (see msh_2d)
%   coeff:      function handle to compute the source function
%   coeff_sp:   function space of the source function
%   coeff_dofs: weitghts of the dofs in the source function
%
% OUTPUT:
%
%   rhs: assembled right-hand side
% 
% Copyright (C) 2011 Rafael Vazquez
% Copyright (C) 2011, 2013 Carlo de Falco
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

function rhs = op_f_v_tp (space, msh, varargin)

  if (numel (varargin) < 2)
      fun_coef = true;
      coeff = varargin{1};
  elseif (numel (varargin) == 2)
      fun_coef = false;      
      [coeff_sp, coeff_dofs] = deal (varargin{:});
  else
      error ('op_f_v_tp: wrong number of input parameters');
  end

  rhs = zeros (space.ndof, 1);

  ndim = numel (msh.qn);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col);

    if (fun_coef)
        for idim = 1:ndim
            x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
        end
        rhs = rhs + op_f_v (sp_col, msh_col, coeff (x{:}));
    else
        rhs = rhs + op_f_v (sp_col, msh_col, ...
                            sp_eval (coeff_dofs, coeff_sp, msh_col));
    end
  end

end
