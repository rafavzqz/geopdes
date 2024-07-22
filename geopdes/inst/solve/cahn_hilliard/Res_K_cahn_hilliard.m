% RES_K_CAHN_HILLIARD: Cahn-Hilliard equation, computation of the residual
% for Newton's method, and of some necessary matrices for the method.
%  It is called from generalized_alpha_step_cahn_hilliard.
%
% Copyright (C) 2023, 2024 Michele Torre, Rafael Vazquez
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

function [Res_gl, stiff_mat] = Res_K_cahn_hilliard (space, msh, mass_mat, lapl_mat, ...
                               bnd_mat, Pen, pen_rhs, u_a, udot_a, mu, dmu)

  % Double well (matrices)
  if (isa (space, 'sp_scalar'))
    [term2, term2K] = op_gradfu_gradv_tp (space, msh, u_a, mu, dmu);
  elseif (isa (space, 'sp_multipatch_C1'))
    [term2, term2K] = op_gradfu_gradv_mp (space, msh, u_a, mu, dmu);
  else
    error ('Not implemented for this kind of space')
  end

  % Residual
  Res_gl = mass_mat*udot_a + term2*u_a  + lapl_mat*u_a;

  % Tangent stiffness matrix (mass is not considered here)
  stiff_mat = term2 + term2K + lapl_mat;

  % In case of Neumann BC, add boundary terms
  if (~isempty(bnd_mat))
    Res_gl = Res_gl - (bnd_mat + bnd_mat.') * u_a + Pen*u_a - pen_rhs;
    stiff_mat = stiff_mat - (bnd_mat + bnd_mat.') + Pen;
  end
end
