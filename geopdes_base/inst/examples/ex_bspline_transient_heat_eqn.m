% EX_BSPLINE_TRANSIENT_HEAT_EQN: Solve a 2d parabolic problem with a B-spline discretization. 
%
% Example to solve the diffusion problem
%
%    du/dt - div ( epsilon(x) grad (u)) = f    in Omega = F((0,1)^2)
%                      epsilon(x) du/dn = g    on Gamma_N
%                                     u = h    on Gamma_D
%
% In order to make it work, a script file with the data must also be executed.
% The available data files for this problem are:
%
% - test_square_heat
% - test_square_heat_exp
%
% Copyright (C) 2009, 2010 Carlo de Falco, Rafael Vazquez
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

% Construct geometry structure
geometry  = geo_load (geo_name);

[knots, zeta] = kntrefine (geometry.nurbs.knots, n_sub, degree, regularity);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_2d_tensor_product (zeta, qn, qw);
msh      = msh_push_forward_2d (msh, geometry);

% Construct space structure
sp       = sp_bspline_2d_phys (knots, degree, msh);

% Precompute the coefficients
x = squeeze (msh.geo_map(1,:,:));
y = squeeze (msh.geo_map(2,:,:));
rho     = reshape (c_mass (x, y), msh.nqn, msh.nel);
epsilon = reshape (c_diff (x, y), msh.nqn, msh.nel);

% Assemble the matrices
mass_mat  = op_u_v (sp, sp, msh, rho);
stiff_mat = op_gradu_gradv (sp, sp, msh, epsilon);

% Find degrees of freedom to be assigned for Dirichlet boundary conditions
u = zeros (sp.ndof, 1);
drchlt_dofs = unique ([sp.boundary(drchlt_sides).dofs]);
int_dofs    = setdiff (1:sp.ndof, drchlt_dofs);

% Compute auxiliary functions to call DASPK
rhs_fun = @(t) op_f_v (sp, msh, reshape (f(x, y, t), msh.nqn, msh.nel));
rhs_nmnn = @(t) sp_nmnn (sp, msh, @(x, y, ind) g(x, y, t, ind), nmnn_sides);

drch = @(t) (sp_drchlt_l2_proj(sp, msh, @(x, y, ind) h(x, y, t, ind), drchlt_sides));

rhs_daspk = @(t) rhs_fun(t)(int_dofs) + rhs_nmnn(t)(int_dofs) - ...
   (stiff_mat(int_dofs, drchlt_dofs) + mass_mat(int_dofs, drchlt_dofs))*drch(t);


FCN = @(x, xdot, t) mass_mat(int_dofs, int_dofs)*xdot + ...
                    stiff_mat(int_dofs, int_dofs)*x - rhs_daspk(t);

% Computation of initial time solution
if (exist ('uex', 'var'))
  u_zero = mass_mat \ op_f_v (sp, msh, reshape (uex(x, y, time_save(1)), msh.nqn, msh.nel));
  udot_zero = mass_mat \ (rhs_fun (time_save(1)) - stiff_mat*u_zero);
else
  u_zero = zeros (sp.ndof, 1);
  udot_zero = zeros (sp.ndof, 1);
end

% DASPK: solution of the transient problem
[y, ydot, istate] = daspk (FCN, u_zero(int_dofs), udot_zero(int_dofs), time_save);

% Postprocessing
error_l2 = zeros(size(time_save));
error_h1 = zeros(size(time_save));
for ii = 1:numel (time_save)
  u(int_dofs) = y(ii, :);
  u(drchlt_dofs) = drch (time_save(ii));
  [eu, F] = sp_eval_2d (u, sp, geometry, vtk_pts);
  [X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
  if (exist ('uex', 'var'))
    uex_time = @(x, y) uex(x, y, time_save(ii));
    err_l2  = sp_l2_error (sp, msh, u, uex_time);
    norm_l2 = sp_l2_error (sp, msh, zeros(sp.ndof,1), uex_time);
    error_l2(ii) = err_l2 / norm_l2;
    if (exist ('graduex', 'var'))
      graduex_time = @(x, y) graduex(x, y, time_save(ii));
      err_h1 = sp_h1_error (sp, msh, u, uex_time, graduex_time);
      norm_h1 = sp_h1_error (sp, msh, zeros(sp.ndof,1), uex_time, graduex_time);
      error_h1(ii) = err_h1 / norm_h1;
    end

    subplot (1,2,2)
    surf (X, Y, uex_time (X,Y))
    title ('exact solution'), axis tight
    subplot(1,2,1)
  end

  surf (X, Y, eu)
  title ('numerical solution'), axis tight
end

%!demo
%! test_square_heat
%! ex_bspline_transient_heat_eqn

%!demo
%! test_square_heat_exp
%! ex_bspline_transient_heat_eqn
