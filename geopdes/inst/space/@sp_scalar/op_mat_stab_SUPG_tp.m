% OP_MAT_STAB_SUPG_tp: assemble the stabilization matrix A = [a(i,j)], exploiting the tensor product structure.
%   a(i,j) =  tau_h * { - ( mu * lapl(u_j) , vel \cdot grad v_i )
%			- ( grad (mu) \cdot grad (u_j), vel \cdot grad (v_i) )
%			+ ( vel \cdot grad u_j, vel \cdot grad v_i ) }
%
%   mat = op_mat_stab_SUPG_tp (spu, spv, msh, mu, grad_mu, vel);
%   [rows, cols, values] = op_mat_stab_SUPG_tp (spu, spv, msh, mu, grad_mu, vel)
%
% INPUT:
%   spu:     object representing the space of trial functions (see sp_scalar)
%   spv:     object representing the space of test functions (see sp_scalar)
%   msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%   mu:      function handle for the diffusion coefficient
%   grad_mu: function handle for the gradient of the diffusion coefficient
%   vel:   function handle for the advection coefficient( vectorial function )
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2011, Carlo de Falco
% Copyright (C) 2011, 2014, Rafael Vazquez
% Copyright (C) 2013, Anna Tagliabue
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

function varargout = op_mat_stab_SUPG_tp (space1, space2, msh, coeff, grad_coeff, vel)
  A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);

    sp1_col = sp_evaluate_col (space1, msh_col, 'value', false, 'gradient', true, 'hessian', true);
    sp2_col = sp_evaluate_col (space2, msh_col, 'value', false, 'gradient', true);

    for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end

    A = A + op_mat_stab_SUPG (sp1_col, sp2_col, msh_col, ...
             coeff (x{:}), grad_coeff (x{:}), vel(x{:}));
  end
   
  if (nargout == 1)
    varargout{1} = A;
  elseif (nargout == 3)
    [rows, cols, vals] = find (A);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end

end
