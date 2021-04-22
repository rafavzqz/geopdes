% OP_KL_BENDING_STRESS: assemble the Kirchhoff-Love bending resultants.
%
%   bending = op_KL_bending_stress (spu, spv, msh, E_coeff, nu_coeff, t_coeff);
%   [bending, Kappa] = op_KL_bending_stress (spu, spv, msh, E_coeff, nu_coeff, t_coeff);
%   [bending, Kappa, T] = op_KL_bending_stress (spu, spv, msh, E_coeff, nu_coeff, t_coeff);
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
%  bending:    bending stress resultant
%  Kappa:      bending strain resultant
%  T:    transformation matrix, useful for stress post-processing
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

function varargout = op_KL_bending_stress (ref_sp_u, ref_sp_v, msh, E_coeff, nu_coeff, t_coeff)

  if (nargin ~= 6)
    error ('compute_KL_bending_stress_tensor: invalid number of input arguments')
  end

  nsh_u = ref_sp_u.nsh_max;
  nsh_v = ref_sp_v.nsh_max;
  nqn = msh.nqn;
  nel = msh.nel;

  A1 = reshape(msh.geo_map_jac(:, 1, :, :), [3, nqn, nel]);
  A2 = reshape(msh.geo_map_jac(:, 2, :, :), [3, nqn, nel]);
  A3 = msh.normal;
  A1xA2 = cross(A1, A2);
  a = permute(repmat(reshape(sqrt(sum(A1xA2 .^ 2)), [nqn, nel]), [1, 1, 3]), [3, 1, 2]);
    
  invJacT = geopdes_invT__(msh.geo_map_jac);
  Act1 = reshape(invJacT(:, 1, :, :), [3, nqn, nel]);
  Act2 = reshape(invJacT(:, 2, :, :), [3, nqn, nel]);

  A11 = reshape(msh.geo_map_der2(:, 1, 1, :, :), [3, nqn, nel]);
  A12 = reshape(msh.geo_map_der2(:, 1, 2, :, :), [3, nqn, nel]);
  A22 = reshape(msh.geo_map_der2(:, 2, 2, :, :), [3, nqn, nel]);

  A11_A3 = permute(reshape(repmat(reshape(sum(A11 .* A3), [nqn, nel]), [3, 1, 1]), [nqn, 3, nel]), [2, 1, 3]);
  A12_A3 = permute(reshape(repmat(reshape(sum(A12 .* A3), [nqn, nel]), [3, 1, 1]), [nqn, 3, nel]), [2, 1, 3]);
  A22_A3 = permute(reshape(repmat(reshape(sum(A22 .* A3), [nqn, nel]), [3, 1, 1]), [nqn, 3, nel]), [2, 1, 3]);

  coef_11 = (A11 - ( A11_A3 .* A3 )) ./ a;
  coef_12 = (A12 - A12_A3 .* A3) ./ a;
  coef_22 = (A22 - A22_A3 .* A3) ./ a;

  coef_1_11 = cross(coef_11, A2);
  coef_1_12 = cross(coef_12, A2);
  coef_1_22 = cross(coef_22, A2);

  coef_2_11 = cross(A1, coef_11);
  coef_2_12 = cross(A1, coef_12);
  coef_2_22 = cross(A1, coef_22);

  A3_v = permute(repmat(A3, [1, 1, 1, nsh_v]), [1, 2, 4, 3]);

  coef_1_11_v = permute(repmat(coef_1_11, [1, 1, 1, nsh_v]), [1, 2, 4, 3]);
  coef_1_12_v = permute(repmat(coef_1_12, [1, 1, 1, nsh_v]), [1, 2, 4, 3]);
  coef_1_22_v = permute(repmat(coef_1_22, [1, 1, 1, nsh_v]), [1, 2, 4, 3]);

  coef_2_11_v = permute(repmat(coef_2_11, [1, 1, 1, nsh_v]), [1, 2, 4, 3]);
  coef_2_12_v = permute(repmat(coef_2_12, [1, 1, 1, nsh_v]), [1, 2, 4, 3]);
  coef_2_22_v = permute(repmat(coef_2_22, [1, 1, 1, nsh_v]), [1, 2, 4, 3]);

  Kappa = zeros(2, 2, nqn, nsh_v, nel);
  Kappa(1, 1, :, :, :) = sum(-A3_v .* reshape(ref_sp_v.shape_function_hessians(:, 1, 1, :, :, :), [3, nqn, nsh_v, nel])) ...
  + sum(coef_1_11_v .* reshape(ref_sp_v.shape_function_gradients(:, 1, :, :, :), [3, nqn, nsh_v, nel])) ...
  + sum(coef_2_11_v .* reshape(ref_sp_v.shape_function_gradients(:, 2, :, :, :), [3, nqn, nsh_v, nel]));

  Kappa(2, 2, :, :, :) = sum(-A3_v .* reshape(ref_sp_v.shape_function_hessians(:, 2, 2, :, :, :), [3, nqn, nsh_v, nel])) ...
  + sum(coef_1_22_v .* reshape(ref_sp_v.shape_function_gradients(:, 1, :, :, :), [3, nqn, nsh_v, nel])) ...
  + sum(coef_2_22_v .* reshape(ref_sp_v.shape_function_gradients(:, 2, :, :, :), [3, nqn, nsh_v, nel]));

  Kappa(1, 2, :, :, :) = sum(-A3_v .* reshape(ref_sp_v.shape_function_hessians(:, 1, 2, :, :, :), [3, nqn, nsh_v, nel])) ...
  + sum(coef_1_12_v .* reshape(ref_sp_v.shape_function_gradients(:, 1, :, :, :), [3, nqn, nsh_v, nel])) ...
  + sum(coef_2_12_v .* reshape(ref_sp_v.shape_function_gradients(:, 2, :, :, :), [3, nqn, nsh_v, nel]));

  Kappa(2, 1, :, :, :) = Kappa(1, 2, :, :, :);
  
  clear A3_v coef_1_11_v coef_1_12_v coef_1_22_v coef_2_11_v coef_2_12_v coef_2_22_v
  
  A3_u = permute(repmat(A3, [1, 1, 1, nsh_u]), [1, 2, 4, 3]);

  coef_1_11_u = permute(repmat(coef_1_11, [1, 1, 1, nsh_u]), [1, 2, 4, 3]);
  coef_1_12_u = permute(repmat(coef_1_12, [1, 1, 1, nsh_u]), [1, 2, 4, 3]);
  coef_1_22_u = permute(repmat(coef_1_22, [1, 1, 1, nsh_u]), [1, 2, 4, 3]);

  coef_2_11_u = permute(repmat(coef_2_11, [1, 1, 1, nsh_u]), [1, 2, 4, 3]);
  coef_2_12_u = permute(repmat(coef_2_12, [1, 1, 1, nsh_u]), [1, 2, 4, 3]);
  coef_2_22_u = permute(repmat(coef_2_22, [1, 1, 1, nsh_u]), [1, 2, 4, 3]);

  aux = zeros(2, 2, nqn, nsh_u, nel);
  aux(1, 1, :, :, :) = sum(-A3_u .* reshape(ref_sp_u.shape_function_hessians(:, 1, 1, :, :, :), [3, nqn, nsh_u, nel])) ...
  + sum(coef_1_11_u .* reshape(ref_sp_u.shape_function_gradients(:, 1, :, :, :), [3, nqn, nsh_u, nel])) ...
  + sum(coef_2_11_u .* reshape(ref_sp_u.shape_function_gradients(:, 2, :, :, :), [3, nqn, nsh_u, nel]));

  aux(2, 2, :, :, :) = sum(-A3_u .* reshape(ref_sp_u.shape_function_hessians(:, 2, 2, :, :, :), [3, nqn, nsh_u, nel])) ...
  + sum(coef_1_22_u .* reshape(ref_sp_u.shape_function_gradients(:, 1, :, :, :), [3, nqn, nsh_u, nel])) ...
  + sum(coef_2_22_u .* reshape(ref_sp_u.shape_function_gradients(:, 2, :, :, :), [3, nqn, nsh_u, nel]));

  aux(1, 2, :, :, :) = sum(-A3_u .* reshape(ref_sp_u.shape_function_hessians(:, 1, 2, :, :, :), [3, nqn, nsh_u, nel])) ...
  + sum(coef_1_12_u .* reshape(ref_sp_u.shape_function_gradients(:, 1, :, :, :), [3, nqn, nsh_u, nel])) ...
  + sum(coef_2_12_u .* reshape(ref_sp_u.shape_function_gradients(:, 2, :, :, :), [3, nqn, nsh_u, nel]));

  aux(2, 1, :, :, :) = aux(1, 2, :, :, :);

  Aco1 = reshape(msh.geo_map_jac(:, 1, :, :), [3, nqn, nel]);

  E1 = reshape(Aco1 ./ sqrt(sum(Aco1 .^ 2)),[3, nqn, nel]);
  E2 = reshape(Act2 ./ sqrt(sum(Act2 .^ 2)),[3, nqn, nel]);
  
  eg11 = reshape(sum(E1 .* Act1), [nqn, nel]);
  eg22 = reshape(sum(E2 .* Act2), [nqn, nel]);
  eg12 = reshape(sum(E1 .* Act2), [nqn, nel]);
  eg21 = reshape(sum(E2 .* Act1), [nqn, nel]);

  T = zeros([3, 3, nqn, nel]);

  T(1,1,:,:) = (eg11 .^ 2);
  T(1,2,:,:) = eg12 .^ 2;
  T(1,3,:,:) = 2 .* eg11 .* eg12;
  T(2,1,:,:) = eg21 .* eg21;
  T(2,2,:,:) = eg22 .* eg22;
  T(2,3,:,:) = 2 .* eg21 .* eg22;
  T(3,1,:,:) = 2 .* eg11 .* eg21;
  T(3,2,:,:) = 2 .* eg12 .* eg22;
  T(3,3,:,:) = 2 .*(eg11 .* eg22 + eg12 .* eg21);
  
  T_transp = permute(T,[2 1 3 4]);
  
  D_bending = zeros([3, 3, nqn, nel]);
  
  coeff = E_coeff .* t_coeff.^3 ./( 12 * (1 - nu_coeff.^2) );
  
  D_bending(1,1,:,:) = coeff;
  D_bending(1,2,:,:) = coeff .* nu_coeff;
  D_bending(2,1,:,:) = coeff .* nu_coeff;
  D_bending(2,2,:,:) = coeff;
  D_bending(3,3,:,:) = coeff .* (1 - nu_coeff)./2;
   
  % this should be vectorized 
  for iqn = 1:nqn
      for iel = 1:nel
        D_bending(:,:,iqn,iel) = T_transp(:,:,iqn,iel) * D_bending(:,:,iqn,iel) * T(:,:,iqn,iel);     
      end
  end

  T_1_1 = permute(repmat(reshape(D_bending(1,1,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);
  T_1_2 = permute(repmat(reshape(D_bending(1,2,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);
  T_1_3 = permute(repmat(reshape(D_bending(1,3,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);
  T_2_1 = permute(repmat(reshape(D_bending(2,1,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);
  T_2_2 = permute(repmat(reshape(D_bending(2,2,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);
  T_2_3 = permute(repmat(reshape(D_bending(2,3,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);
  T_3_1 = permute(repmat(reshape(D_bending(3,1,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);
  T_3_2 = permute(repmat(reshape(D_bending(3,2,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);
  T_3_3 = permute(repmat(reshape(D_bending(3,3,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);

  bending = zeros(2, 2, nqn, nsh_u, nel);
    
  bending(1, 1, :, :, :) = T_1_1 .* reshape(aux(1, 1, :, :, :), [nqn, nsh_u, nel]) + ...
            T_1_2 .* reshape(aux(2, 2, :, :, :), [nqn, nsh_u, nel]) + ...
            T_1_3 .* reshape(aux(1, 2, :, :, :), [nqn, nsh_u, nel]);

  bending(1, 2, :, :, :) = .5 * (T_3_1 .* reshape(aux(1, 1, :, :, :), [nqn, nsh_u, nel]) + ...
            T_3_2 .* reshape(aux(2, 2, :, :, :), [nqn, nsh_u, nel]) + ...
            T_3_3 .* reshape(aux(1, 2, :, :, :), [nqn, nsh_u, nel]));

  bending(2, 2, :, :, :) = T_2_1 .* reshape(aux(1, 1, :, :, :), [nqn, nsh_u, nel]) + ...
            T_2_2 .* reshape(aux(2, 2, :, :, :), [nqn, nsh_u, nel]) + ...
            T_2_3 .* reshape(aux(1, 2, :, :, :), [nqn, nsh_u, nel]);

  bending(2, 1, :, :, :) =   bending(1, 2, :, :, :);

  clear aux
  
  if (nargout == 1)
    varargout{1} = bending;
  elseif (nargout == 2)
    varargout{1} = bending;
    varargout{2} = Kappa;
  else
    varargout{1} = bending;
    varargout{2} = Kappa;
    varargout{3} = T;
  end

end

