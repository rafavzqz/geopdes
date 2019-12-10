% OP_KL_SHELLS: assemble the Kirchhoff-Love stiffness matrix K.
%
%   mat = op_KL_shells (spu, spv, msh, E_coeff, nu_coeff, t_coeff);
%
% INPUT:
%
%  spu:   structure representing the space of trial functions (see sp_vector/sp_evaluate_col)
%  spv:   structure representing the space of test functions  (see sp_vector/sp_evaluate_col)
%  msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%  E_coeff:  coefficients for the Young's modulus
%  nu_coeff: coefficients for the Poisson's ratio
%  t_coeff:  thickness of the shell
%
% OUTPUT:
%
%  mat:    assembled stiffness matrix
% 
% Copyright (C) 2018, 2019 Pablo Antolin, Luca Coradello
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

function mat = op_KL_shells (ref_sp_u, ref_sp_v, msh, E_coeff, nu_coeff, t_coeff)

   wq = msh.jacdet .* msh.quad_weights;

   [bending_stress, Kappa] = op_KL_bending_stress (ref_sp_u, ref_sp_v, msh, E_coeff .* wq, nu_coeff, t_coeff);
   [membrane_stress, Epsilon] = op_KL_membrane_stress (ref_sp_u, ref_sp_v, msh, E_coeff .* wq, nu_coeff, t_coeff);
   
   bending_stress = reshape(bending_stress, [], msh.nqn, ref_sp_u.nsh_max, msh.nel);
   membrane_stress = reshape(membrane_stress, [], msh.nqn, ref_sp_u.nsh_max, msh.nel);

   Kappa = reshape(Kappa, [], msh.nqn, ref_sp_v.nsh_max, msh.nel);
   Epsilon = reshape(Epsilon, [], msh.nqn, ref_sp_v.nsh_max, msh.nel);

   rows = zeros (msh.nel * ref_sp_u.nsh_max * ref_sp_v.nsh_max, 1);
   cols = zeros (msh.nel * ref_sp_u.nsh_max * ref_sp_v.nsh_max, 1);
   values = zeros (msh.nel * ref_sp_u.nsh_max * ref_sp_v.nsh_max, 1);
  
   ncounter = 0;
   
   for iel = 1:msh.nel
     if (all (msh.jacdet(:,iel)))

       ndof_u = ref_sp_u.nsh_max;
       ndof_v = ref_sp_v.nsh_max;

       bending_stress_iel = reshape(bending_stress (:, :, :, iel), [], msh.nqn, 1, ndof_u);
       membrane_stress_iel = reshape(membrane_stress (:, :, :, iel), [], msh.nqn, 1, ndof_u);
       
       Kappa_iel = reshape(Kappa (:, :, :, iel), [], msh.nqn, ndof_v, 1);
       Epsilon_iel = reshape(Epsilon (:, :, :, iel), [], msh.nqn, ndof_v, 1);
       
       tmp1 = sum (bsxfun (@times, bending_stress_iel, Kappa_iel), 1);
       tmp2 = sum (bsxfun (@times, membrane_stress_iel, Epsilon_iel), 1);
       
       elementary_values =   reshape (sum (tmp1, 2), ndof_v, ndof_u) + ...
                            reshape (sum (tmp2, 2), ndof_v, ndof_u);
       [rows_loc, cols_loc] = ndgrid (ref_sp_v.connectivity(:,iel), ref_sp_u.connectivity(:,iel));
       indices = rows_loc & cols_loc;
       rows(ncounter+(1:ref_sp_u.nsh(iel)*ref_sp_v.nsh(iel))) = rows_loc(indices);
       cols(ncounter+(1:ref_sp_u.nsh(iel)*ref_sp_v.nsh(iel))) = cols_loc(indices);
       values(ncounter+(1:ref_sp_u.nsh(iel)*ref_sp_v.nsh(iel))) = elementary_values(indices);
       ncounter = ncounter + ref_sp_u.nsh(iel)*ref_sp_v.nsh(iel);

%%%%%%%%%%%%%%%% old version, works when test and trial space are the same
%        % Here, the loops over the two sets of functions and the loop over the quadrature points are
%        % vectorized into a single operation.
%        ndof = ref_sp_u.nsh_max;
%        mat(ref_sp_v.connectivity(:, iel), ref_sp_u.connectivity(:, iel)) = ...
%            mat(ref_sp_v.connectivity(:, iel), ref_sp_u.connectivity(:, iel)) ...
%            + reshape(sum(sum(sum(repmat(bending_stress (:, :, :, :, iel), [1, 1, 1, 1, ndof]) ...
%                            .* permute(repmat(Kappa (:, :, :, :, iel), [1, 1, 1, 1, ndof]), [1, 2, 3, 5, 4])))) ...
%                  + sum(sum(sum(repmat(membrane_stress (:, :, :, :, iel), [1, 1, 1, 1, ndof]) ...
%                            .* permute(repmat(Epsilon (:, :, :, :, iel), [1, 1, 1, 1, ndof]), [1, 2, 3, 5, 4])))), [ndof, ndof]);

     else
       warning ('geopdes:jacdet_zero_at_quad_node',...
           'op_KL: singular map in element number %d', iel)
     end
   end
   mat = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), ref_sp_v.ndof, ref_sp_u.ndof);
end

