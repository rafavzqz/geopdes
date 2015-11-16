% OP_RHS_STAB_SUPG_TP: assemble the rhs stabilization vector rhs_auxiliar = [r(i)],
% r(i) = tau_h * { ( f, vel \cdot grad v_i ) }, exploiting the tensor product structure.
%
%   rhs = op_rhs_stab_SUPG_tp (space, msh, mu, vel, f);
%
% INPUT:
%
%   spv:   object representing the function space (see sp_scalar)
%   msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%   mu:    function handle for the diffusion coefficient
%   vel:   function handle for the advection coefficient( vectorial function )
%   f:     function handle to compute the source function
%
% OUTPUT:
%
%   rhs_auxiliar: assembled right-hand side relative to the stabilization terms to be added to rhs
% 
% Copyright (C) 2011, Carlo de Falco
% Copyright (C) 2011, 2014 Rafael Vazquez
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

function rhs_auxiliar = op_rhs_stab_SUPG_tp (space, msh, coeff_mu, vel, f)

  rhs_auxiliar = zeros (space.ndof, 1);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col = sp_evaluate_col (space, msh_col, 'value', false, 'gradient', true);

    for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end

    rhs_auxiliar = rhs_auxiliar + op_rhs_stab_SUPG (sp_col, msh_col, coeff_mu (x{:}), vel (x{:}), f (x{:}));

  end

end
