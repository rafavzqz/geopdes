% OP_KL_MEMBRANE_STRESS: assemble the Kirchhoff-Love membrane
% resultants
%
%   membrane = op_KL_membrane_stress (spu, spv, msh, E_coeff, nu_coeff, t_coeff);
%   [membrane, Epsilon] = op_KL_membrane_stress (spu, spv, msh, E_coeff, nu_coeff, t_coeff);
%   [membrane, Epsilon T] = op_KL_membrane_stress (spu, spv, msh, E_coeff, nu_coeff, t_coeff);
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
%  membrane: membrane stress resultant
%  Epsilon:  membrane strain resultant
%  T:        transformation matrix, useful for stress post-processing
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

function varargout = op_KL_membrane_stress (ref_sp_u, ref_sp_v, msh, E_coeff, nu_coeff, t_coeff)

  if (nargin ~= 6)
    error ('compute_KL_membrane_stress_tensor: invalid number of input arguments')
  end

  nsh_u = ref_sp_u.nsh_max;
  nsh_v = ref_sp_v.nsh_max;
  nqn = msh.nqn;
  nel = msh.nel;

  A1_v = permute(reshape(repmat(reshape(msh.geo_map_jac(:, 1, :, :), [3, nqn, nel]), [1, 1, 1, nsh_v]), [3, msh.nqn, msh.nel, nsh_v]), [1, 2, 4, 3]);
  A2_v = permute(reshape(repmat(reshape(msh.geo_map_jac(:, 2, :, :), [3, nqn, nel]), [1, 1, 1, nsh_v]), [3, msh.nqn, msh.nel, nsh_v]), [1, 2, 4, 3]);

  invJacT = geopdes_invT__(msh.geo_map_jac);
  Act1 = reshape(invJacT(:, 1, :, :), [3, nqn, nel]);
  Act2 = reshape(invJacT(:, 2, :, :), [3, nqn, nel]);

  g1_v = reshape(ref_sp_v.shape_function_gradients(:, 1, :, :, :), [3, nqn, nsh_v, nel]);
  g2_v = reshape(ref_sp_v.shape_function_gradients(:, 2, :, :, :), [3, nqn, nsh_v, nel]);

  Epsilon = zeros(2, 2, msh.nqn, nsh_v, msh.nel);
  Epsilon(1, 1, :, :, :) = sum(A1_v .* g1_v);
  Epsilon(2, 2, :, :, :) = sum(A2_v .* g2_v);
  Epsilon(1, 2, :, :, :) = 0.5 * (sum(A2_v .* g1_v) + sum(A1_v .* g2_v));
  Epsilon(2, 1, :, :, :) = Epsilon(1, 2, :, :, :);

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
  
  D_membrane = zeros([3, 3, nqn, nel]);
  
  coeff = E_coeff .* t_coeff ./ (1 - nu_coeff.^2);
  
  D_membrane(1,1,:,:) = coeff;
  D_membrane(1,2,:,:) = coeff .* nu_coeff;
  D_membrane(2,1,:,:) = coeff .* nu_coeff;
  D_membrane(2,2,:,:) = coeff;
  D_membrane(3,3,:,:) = coeff .* (1 - nu_coeff)./2;
   
  % this should be vectorized 
  for iqn = 1:nqn
      for iel = 1:nel
        D_membrane(:,:,iqn,iel) = T_transp(:,:,iqn,iel) * D_membrane(:,:,iqn,iel) * T(:,:,iqn,iel);     
      end
  end

  T_1_1 = permute(repmat(reshape(D_membrane(1,1,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);
  T_1_2 = permute(repmat(reshape(D_membrane(1,2,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);
  T_1_3 = permute(repmat(reshape(D_membrane(1,3,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);
  T_2_1 = permute(repmat(reshape(D_membrane(2,1,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);
  T_2_2 = permute(repmat(reshape(D_membrane(2,2,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);
  T_2_3 = permute(repmat(reshape(D_membrane(2,3,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);
  T_3_1 = permute(repmat(reshape(D_membrane(3,1,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);
  T_3_2 = permute(repmat(reshape(D_membrane(3,2,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);
  T_3_3 = permute(repmat(reshape(D_membrane(3,3,:,:), [nqn, nel]), [1, 1, nsh_u]), [1, 3, 2]);

  A1_u = permute(reshape(repmat(reshape(msh.geo_map_jac(:, 1, :, :), [3, nqn, nel]), [1, 1, 1, nsh_u]), [3, msh.nqn, msh.nel, nsh_u]), [1, 2, 4, 3]);
  A2_u = permute(reshape(repmat(reshape(msh.geo_map_jac(:, 2, :, :), [3, nqn, nel]), [1, 1, 1, nsh_u]), [3, msh.nqn, msh.nel, nsh_u]), [1, 2, 4, 3]);

  g1_u = reshape(ref_sp_u.shape_function_gradients(:, 1, :, :, :), [3, nqn, nsh_u, nel]);
  g2_u = reshape(ref_sp_u.shape_function_gradients(:, 2, :, :, :), [3, nqn, nsh_u, nel]);

  aux = zeros(2, 2, msh.nqn, nsh_u, msh.nel);
  aux(1, 1, :, :, :) = sum(A1_u .* g1_u);
  aux(2, 2, :, :, :) = sum(A2_u .* g2_u);
  aux(1, 2, :, :, :) = 0.5 * (sum(A2_u .* g1_u) + sum(A1_u .* g2_u));
  aux(2, 1, :, :, :) = aux(1, 2, :, :, :);
  
  membrane = zeros(2, 2, nqn, nsh_u, nel);
    
  membrane(1, 1, :, :, :) = T_1_1 .* reshape(aux(1, 1, :, :, :), [nqn, nsh_u, nel]) + ...
            T_1_2 .* reshape(aux(2, 2, :, :, :), [nqn, nsh_u, nel]) + ...
            T_1_3 .* reshape(aux(1, 2, :, :, :), [nqn, nsh_u, nel]);

  membrane(1, 2, :, :, :) = .5 * (T_3_1 .* reshape(aux(1, 1, :, :, :), [nqn, nsh_u, nel]) + ...
            T_3_2 .* reshape(aux(2, 2, :, :, :), [nqn, nsh_u, nel]) + ...
            T_3_3 .* reshape(aux(1, 2, :, :, :), [nqn, nsh_u, nel]));

  membrane(2, 2, :, :, :) = T_2_1 .* reshape(aux(1, 1, :, :, :), [nqn, nsh_u, nel]) + ...
            T_2_2 .* reshape(aux(2, 2, :, :, :), [nqn, nsh_u, nel]) + ...
            T_2_3 .* reshape(aux(1, 2, :, :, :), [nqn, nsh_u, nel]);

  membrane(2, 1, :, :, :) =   membrane(1, 2, :, :, :);
  
  clear aux
  
  if (nargout == 1)
    varargout{1} = membrane;
  elseif (nargout == 2)
    varargout{1} = membrane;
    varargout{2} = Epsilon;
  else
    varargout{1} = membrane;
    varargout{2} = Epsilon;
    varargout{3} = T;
  end
  
end

