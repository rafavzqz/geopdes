% OP_KL_SHELLS_MP: assemble the Kirchhoff-Love shell stiffness matrix.
%
%   mat = op_KL_shells_mp (spu, spv, msh, E_coeff, nu_coeff, t_coeff, [patches]);
%   [rows, cols, values] = op_KL_shells_mp (spu, spv, msh, E_coeff, nu_coeff, t_coeff, [patches]);
%
% INPUT:
%
%  spu:   object representing the space of trial functions (see sp_multipatch_C1)
%  spv:   object representing the space of test functions (see sp_multipatch_C1)
%  msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%  E_coeff:  function handle to compute the Young's modulus
%  nu_coeff: function handle to compute the Poisson's ratio
%  t_coeff:  thickness of the shell, scalar value
%  patches:  list of patches where the integrals have to be computed. By default, all patches are selected.

% OUTPUT:
%
%  mat:    assembled stiffness matrix
% 
% Copyright (C) 2022-2023 Rafael Vazquez
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

function A = op_KL_shells_mp (space_u, space_v, msh, E_coeff, nu_coeff, t_coeff, patch_list)

  if (nargin < 7)
    patch_list = 1:msh.npatch;
  end

  if ((space_u.npatch ~= space_v.npatch) || (space_u.npatch ~= msh.npatch))
    error ('op_KL_shells_mp: the number of patches does not coincide')
  end

  A = sparse (msh.rdim*space_v.ndof, msh.rdim*space_u.ndof);

  for iptc = patch_list
    Ap = op_KL_shells_tp (space_u.sp_patch{iptc}, space_v.sp_patch{iptc}, msh.msh_patch{iptc}, E_coeff, nu_coeff, t_coeff);
    
    [Cpatch_u, cols_u] = sp_compute_Cpatch_vector (space_u, iptc, msh.rdim);
    [Cpatch_v, cols_v] = sp_compute_Cpatch_vector (space_v, iptc, msh.rdim);

    A(cols_v,cols_u) = A(cols_v,cols_u) + Cpatch_v.' * Ap * Cpatch_u;
  end
  
end
