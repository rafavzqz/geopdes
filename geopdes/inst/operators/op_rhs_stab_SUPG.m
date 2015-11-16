% OP_RHS_STAB_SUPG: assemble the rhs stabilization vector rhs = [r(i)],
%   rhs(i) = tau_h * { ( f, vel \cdot grad v_i ) }
%
% The current version only works for the scalar case.
%
%   rhs = op_rhs_stab_SUPG (space, msh, mu, vel, f);
%
% INPUT:
%
%   space: structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   mu:    diffusion coefficient evaluated at the quadrature points
%   vel: advection coefficient( vectorial function ), evaluated at the quadrature points
%   coeff: source function evaluated at the quadrature points
%
% OUTPUT:
%
%   rhs: assembled right-hand side for stabilization 
% 
% Copyright (C) 2009, 2010 Carlo de Falco
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

function rhs = op_rhs_stab_SUPG (space, msh, coeff_mu, vel, coeff)

  rhs = zeros (space.ndof, 1);

  coeff = reshape (coeff, space.ncomp, msh.nqn, msh.nel);
  coeff_mu = reshape (coeff_mu, 1, msh.nqn, msh.nel);

  gradv = reshape (space.shape_function_gradients, space.ncomp, [], ...
		   msh.nqn, space.nsh_max, msh.nel);

  ndir = size (gradv, 2);

  p = max (space.degree(:));

  for iel = 1:msh.nel  
    if (all (msh.jacdet(:, iel)))
      vel_iel = reshape (vel(:, :, iel), [], msh.nqn);

% compute parameters relative to the stabilization coefficient
      h_iel = msh.element_size(iel);
      max_coeff = max (abs (coeff_mu(1, :, iel)));
      [max_vel, ind] = max (sqrt (sum (vel_iel.^2., 1)));
% Length in the direction of the velocity. This could be improved.
      h_iel = h_iel / max (abs (vel_iel(:,ind)) / max_vel);
      
      Pe = max_vel * h_iel / (2. * max_coeff);
      tau = h_iel / (2. * max_vel) * min (1., Pe / ( 3. * p * p ));
      
      %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      gradv_iel = permute (gradv(:, :, :, 1:space.nsh(iel), iel), [1 2 4 3]);
      gradv_iel = reshape (gradv_iel, space.ncomp * ndir, space.nsh(iel), msh.nqn);
      gradv_iel = permute (gradv_iel, [1 3 2]);

      %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      jacdet_weights = reshape (msh.jacdet(:, iel) .* msh.quad_weights(:, iel) * tau, 1, msh.nqn);
      coeff_times_jw = bsxfun (@times, jacdet_weights, coeff(:,:,iel));
      jacdet_weights_vel = bsxfun (@times, coeff_times_jw, vel_iel);

      aux_val = bsxfun (@times, jacdet_weights_vel, gradv_iel);

      rhs_loc = sum (sum (aux_val, 1), 2);
      rhs(space.connectivity(1:space.nsh(iel), iel)) = rhs(space.connectivity(1:space.nsh(iel), iel)) + rhs_loc(:); 
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_f_v: singular map in element number %d', iel)
    end
  end

end

