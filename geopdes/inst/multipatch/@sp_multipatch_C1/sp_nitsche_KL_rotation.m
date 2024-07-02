% SP_NITSCHE_KL_ROTATION: impose zero rotation for Kirchhoff-Love shells using Nitsche method
%
%   N_mat = sp_nitsche_KL_rotation (space, msh, bnd_sides, E_coeff, nu_coeff, thickness, C_penalty)
%
% INPUT:
%     
%    space:     space object (see sp_multipatch_C1)
%    msh:       mesh object (see msh_multipatch)
%    bnd_sides: boundary sides on which the rotation free condition is imposed
%    E_coeff:   function handle for the Young modulus
%    nu_coeff:  function handle for the Poisson ratio
%    thickness: scalar value for the thickness
%    C_penalty: parameter for the penalization term
%   
% OUTPUT:
%
%     N_mat:       the computed matrix, to be added in the left hand-side
%
% Copyright (C) 2023, 2024 Giuliano Guarino, Rafael Vazquez
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

function A = sp_nitsche_KL_rotation (space, msh, refs, E_coeff, nu_coeff, thickness, penalty_coeff)

  A = spalloc (msh.rdim*space.ndof, msh.rdim*space.ndof, 3*msh.rdim*space.ndof);

% Compute the matrices to impose the tangential boundary condition weakly
  penalty_coeff = penalty_coeff * max (msh.msh_patch{1}.nel_dir);
  for iref = refs
    for bnd_side = 1:msh.boundaries(iref).nsides
      iptc = msh.boundaries(iref).patches(bnd_side);
      iside = msh.boundaries(iref).faces(bnd_side);

      msh_side = msh_eval_boundary_side (msh.msh_patch{iptc}, iside);
      msh_side_from_interior = msh_boundary_side_from_interior (msh.msh_patch{iptc}, iside);
      
      sp_bnd = space.sp_patch{iptc}.constructor (msh_side_from_interior);
      msh_side_fi = msh_precompute (msh_side_from_interior);
      sp_bnd_param = sp_precompute_param (sp_bnd, msh_side_fi, 'value', true, 'gradient', true, 'hessian', true);
      sp_bnd_param = sp_scalar2vector_param (sp_bnd_param, msh_side_fi, 'value', true, 'gradient', true, 'hessian', true);

      for idim = 1:msh.rdim
        x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
      end

% Evaluating parameters on the boundary
      E_bnd  = reshape(E_coeff  (x{:}), msh_side.nqn, msh_side.nel);
      nu_bnd = reshape(nu_coeff (x{:}), msh_side.nqn, msh_side.nel);

      A_side = op_nitsche_KL_boundary (sp_bnd_param, sp_bnd_param, msh_side, msh_side_fi, E_bnd, nu_bnd, thickness, penalty_coeff);

      [Cpatch_u, cols_u] = sp_compute_Cpatch_vector (space, iptc, msh.rdim);
      A(cols_u,cols_u) = A(cols_u,cols_u) + Cpatch_u.' * A_side * Cpatch_u;
      % rhs(cols_u) = rhs(cols_u) + Cpatch_u.' * rhs_side;
    end
  end

end
