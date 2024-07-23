% SOLVE_CAHN_HILLIARD: solve the Cahn-Hilliard equation, with a generalized alpha discretization in time.
%
% The functions solves the problem of finding u such that
%
%  du/dt - Delta (mu(u) - lambda*Delta u) = 0
%
% with Delta the Laplacian, and mu(u) = alpha u^3 - beta u, and periodic or Neumann boundary conditions.
%
% The values of mu and its derivative are given in problem_data, and used in op_gradfu_gradv_tp.
%
% For details on the problem and the formulation, see
%  H. Gomez, V.M. Calo, Y. Bazilevs, T.J.R. Hughes, CMAME 197 (2008), 4333-4352.
%  H. Gomez, A. Reali, G. Sangalli, J. Comput. Physics 262 (2014), 153-171.
%
% USAGE:
%
%   [geometry, msh, space, results] = solve_cahn_hilliard (problem_data, method_data, save_info)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - periodic_directions: parametric directions along which to apply periodic conditions (may be empty)
%    - lambda:       parameter representing the length scale of the problem, and the width of the interface
%    - mu:           function handle to compute mu (from the double well function)
%    - dmu:          function handle to compute the derivative of mu
%    - Time_max:     final time
%    - fun_u:        initial condition. Equal to zero by default.
%    - fun_udot:     initial condition for time derivative. Equal to zero by default.
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%    - dt:         time step size for generalized-alpha method
%    - rho_inf_gen_alpha: parameter in [0,1], which governs numerical damping of the generalized alpha method
%    - Cpen_nitsche: penalization parameter for Nitsche's method, to impose Neumann conditions
%    - Cpen_projection: penalization parameter to impose zero flux for the initial condition
%
%  save_info: time values at which the solution should be saved (see "results" in the output)
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_cartesian)
%  space:    space object that defines the discrete space (see sp_scalar)
%  results:  a struct with the saved results, containing the following fields:
%    - time: (array of length Ntime) time at which the solution was saved
%    - u:    (size ndof x Ntime) degrees of freedom for the solution
%    - udot: (size ndof x Ntime) degrees of freedom for the time derivative
%
% Only periodic and Neumann boundary conditions are implemented. Neumann
%  conditions are considered by default.
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

function [geometry, msh, space, results] = solve_cahn_hilliard (problem_data, method_data, save_info)

%%-------------------------------------------------------------------------
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

%%-------------------------------------------------------------------------
% Construct geometry structure
geometry  = geo_load (geo_name);
[knots, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);

%%-------------------------------------------------------------------------
% Check for periodic conditions, and consistency with other boundary conditions
if (exist('periodic_directions', 'var'))
  knots = kntunclamp (knots, degree, regularity, periodic_directions);
else
  periodic_directions = [];
end

if (exist ('nmnn_sides','var') && ~isempty (nmnn_sides))
  disp('User defined Neumann sides deleted')
  clear nmnn_sides
end

%%-------------------------------------------------------------------------
% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);
  
% Construct space structure
space    = sp_bspline (knots, degree, msh, [], periodic_directions);

%%-------------------------------------------------------------------------
% Generalized alpha parameters
a_m = .5*((3-rho_inf_gen_alpha)/(1+rho_inf_gen_alpha));
a_f =  1/(1+rho_inf_gen_alpha);
gamma =  .5 + a_m - a_f;

%%-------------------------------------------------------------------------
% No flux b.c. (essential boundary condition)
% Set Neumann boundary conditions for non-periodic sides
nmnn_sides   = [];
for idir = 1:msh.ndim
  if (~ismember(idir, periodic_directions))
    nmnn_sides = [nmnn_sides, 2*(idir-1)+[1 2]];
  end
end

%%-------------------------------------------------------------------------
% Precompute matrices

% Compute the mass matrix
mass_mat = op_u_v_tp (space,space,msh);

% Compute the laplace matrix
lapl_mat = op_laplaceu_laplacev_tp (space, space, msh, lambda);

% Compute the boundary term
bnd_mat = op_nitsche_consistency_cahn_hilliard (space, msh, nmnn_sides, lambda);

% Compute the penalty matrix
%[Pen, pen_rhs] = penalty_matrix (space, msh, nmnn_sides, Cpen_nitsche);
[Pen, pen_rhs] = op_penalty_dudn (space, msh, nmnn_sides, Cpen_nitsche);

%%-------------------------------------------------------------------------
% Initial conditions
time = problem_data.initial_time;
if (exist('fun_u', 'var') && ~isempty(fun_u))
  rhs = op_f_v_tp (space, msh, fun_u);
  u_n = (mass_mat + Cpen_projection/Cpen_nitsche * Pen)\rhs;
else
  u_n = zeros(space.ndof, 1);
end

if (exist('fun_udot', 'var') && ~isempty(fun_udot))
  rhs = op_f_v_tp (space, msh, fun_udot);
  udot_n = (mass_mat + Cpen_projection/Cpen_nitsche * Pen)\rhs;
else
  udot_n = zeros(space.ndof, 1);
end

%%-------------------------------------------------------------------------
% Initialize structure to store the results
save_info = save_info(save_info>=problem_data.initial_time & save_info<=problem_data.Time_max);

results.u = zeros(length(u_n), length(save_info));
results.udot = zeros(length(u_n), length(save_info));
results.time = zeros(length(save_info), 1);
save_info(end+1) = problem_data.Time_max + 1e5;

% Save initial conditions
save_id = 1;
if (time >= save_info(1))
  results.u(:,save_id) = u_n;
  results.udot(:,save_id) = udot_n;
  results.time(save_id) = time;
  save_id = 2;
  save_info = save_info(save_info > time);
end

%%-------------------------------------------------------------------------
% Loop over time steps
while time < Time_max
  disp('----------------------------------------------------------------')
  disp(strcat('time step t=',num2str(time)))

  [u_n1, udot_n1] = generalized_alpha_step_cahn_hilliard (u_n, udot_n, dt, a_m, a_f, gamma, mu, dmu, ...
                    mass_mat, lapl_mat, bnd_mat, Pen, pen_rhs, space, msh);

  % Time step update
  time = time + dt;
  u_n = u_n1;
  udot_n = udot_n1;

  % Check max time
  if (time + dt > Time_max)
    dt = Time_max-time;
  end

  % Store results
  if (time >= save_info(1))
    results.u(:,save_id) = u_n;
    results.udot(:,save_id) = udot_n;
    results.time(save_id) = time;
    save_id = save_id + 1;
    save_info = save_info(save_info > time);
  end
end

disp('----------------------------------------------------------------')
disp(strcat('END ANALYSIS t=',num2str(time)))
disp('----------------------------------------------------------------')

% Crop results
results.u = results.u(:,1:save_id-1);
results.udot = results.udot(:,1:save_id-1);
results.time = results.time(1:save_id-1);

end

%--------------------------------------------------------------------------
% FUNCTIONS
%--------------------------------------------------------------------------


% %--------------------------------------------------------------------------
% % One step of generalized alpha-method
% %--------------------------------------------------------------------------
% 
% function [u_n1, udot_n1] = generalized_alpha_step(u_n, udot_n, dt, a_m, a_f, gamma, mu, dmu, ...
%                        mass_mat, lapl_mat, bnd_mat, Pen, pen_rhs, space, msh)
% 
% % Convergence criteria
%   n_max_iter = 20;
%   tol_rel_res = 1e-10;
%   tol_abs_res = 1e-10;
% 
% % Predictor step
%   u_n1 = u_n;
%   udot_n1 = (gamma-1)/gamma * udot_n; 
% 
% % Newton loop
%   for iter = 0:n_max_iter
% 
%   % Field at alpha level
%     udot_a = udot_n + a_m *(udot_n1-udot_n);
%     u_a = u_n + a_f *(u_n1-u_n);
% 
%   % Compute the residual (internal)
%     [Res_gl, stiff_mat] = Res_K_cahn_hilliard (space, msh, ...
%                           mass_mat, lapl_mat, bnd_mat, Pen, pen_rhs, ...
%                           u_a, udot_a, mu, dmu);
% 
%   % Convergence check
%     if (iter == 0)
%       norm_res_0 = norm(Res_gl);
%     end
%     norm_res = norm(Res_gl);
% 
%     if (norm_res/norm_res_0 < tol_rel_res)
%       disp(strcat('iteration n°=',num2str(iter)))
%       disp(strcat('norm (abs) residual=',num2str(norm_res)))
%       break
%     end
%     if (norm_res<tol_abs_res)
%       disp(strcat('iteration n°=',num2str(iter)))
%       disp(strcat('norm absolute residual=',num2str(norm_res)))
%       break
%     end
%     if (iter == n_max_iter)
%       disp(strcat('Newton reached the maximum number of iterations'))
%       disp(strcat('norm residual=',num2str(norm_res)))
%     end
% 
%   % Compute the update, and update the solution
%     A_gl = a_m * mass_mat + a_f * gamma * dt * stiff_mat ; 
%     d_udot = - A_gl\Res_gl;
% 
%     udot_n1 = udot_n1 + d_udot;
%     u_n1 = u_n1 + gamma * dt * d_udot;
%   end
% 
% end
% 
% %--------------------------------------------------------------------------
% % Cahn-Hilliard residual and tangent matrix
% %--------------------------------------------------------------------------
% 
% function [Res_gl, stiff_mat] = Res_K_cahn_hilliard(space, msh, ...
%                           mass_mat, lapl_mat, bnd_mat, Pen, pen_rhs, u_a, udot_a, mu, dmu)
% 
%   % Double well (matrices)
%   [term2, term2K] = op_gradfu_gradv_tp (space, msh, u_a, mu, dmu);
% 
%   % Residual
%   Res_gl = mass_mat*udot_a + term2*u_a  + lapl_mat*u_a;
% 
%   % Tangent stiffness matrix (mass is not considered here)
%   stiff_mat = term2 + term2K + lapl_mat;
% 
%   % In case of Neumann BC, add boundary terms
%   if (~isempty(bnd_mat))
%     Res_gl = Res_gl - (bnd_mat + bnd_mat.') * u_a + Pen*u_a - pen_rhs;
%     stiff_mat = stiff_mat - (bnd_mat + bnd_mat.') + Pen;
%   end
% end
% 
% %--------------------------------------------------------------------------
% % Boundary term, \int_\Gamma (\Delta u) (\partial v / \partial n)
% %--------------------------------------------------------------------------
% 
% function [A] = int_boundary_term (space, msh,  lambda, nmnn_sides)
% 
%   if (~isempty(nmnn_sides))
% 
%     A =  spalloc (space.ndof, space.ndof, 3*space.ndof);
% 
%     for iside = 1:numel(nmnn_sides)   
% 
%       msh_side = msh_eval_boundary_side (msh, nmnn_sides(iside) ) ;
%       msh_side_int = msh_boundary_side_from_interior (msh, nmnn_sides(iside) ) ;
%       sp_side = space.constructor ( msh_side_int) ;
%       sp_side = sp_precompute ( sp_side , msh_side_int , 'gradient' , true, 'laplacian', true );
% 
%       for idim = 1:msh.rdim
%         x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
%       end
%       coe_side = lambda (x{:});
% 
%       A =  A + op_gradv_n_laplaceu(sp_side ,sp_side ,msh_side, coe_side);
%     end
% 
%   else
%     A = [];
%   end
% end

% %--------------------------------------------------------------------------
% % Penalty term
% %--------------------------------------------------------------------------
% 
% function [P, rhs] = penalty_matrix (space, msh, nmnn_sides, pen)
% 
%   P =  spalloc (space.ndof, space.ndof, 3*space.ndof);
%   rhs = zeros(space.ndof,1);
% 
%   for iside = 1:numel(nmnn_sides)
%     [mass_pen, rhs_pen] = penalty_grad (space, msh, nmnn_sides(iside), pen);
%     P = P + mass_pen; 
%     rhs = rhs + rhs_pen;
%   end
% 
% end
% 
% 
% function [mass_pen, rhs_pen] = penalty_grad (space, msh, i_side, pen)
% 
% msh_side = msh_eval_boundary_side (msh, i_side);
% msh_side_int = msh_boundary_side_from_interior (msh, i_side);
% sp_side = space.constructor (msh_side_int);
% sp_side = sp_precompute (sp_side , msh_side_int , 'gradient', true);
% 
% coe_side = pen .* msh_side.charlen; 
% mass_pen = op_gradu_n_gradv_n(sp_side, sp_side, msh_side, coe_side);
% 
% rhs_pen  = zeros(space.ndof,1);
% 
% end

%--------------------------------------------------------------------------
% Check flux through the boundaries
%--------------------------------------------------------------------------

function flux = check_flux_phase_field(space, msh, uhat, uhat0, sides)

flux = 0;

for iside = 1:numel(sides)   

  msh_side = msh_eval_boundary_side (msh, sides(iside));
  msh_side_int = msh_boundary_side_from_interior (msh, sides(iside));
  sp_side = space.constructor (msh_side_int);
  sp_side = sp_precompute (sp_side , msh_side_int , 'gradient', true );

  gradu = sp_eval_msh (uhat-uhat0, sp_side, msh_side, 'gradient');
    
  valu = zeros(sp_side.ncomp, size(msh_side.quad_weights,1), size(msh_side.quad_weights,2));
  for idim = 1:msh.rdim
    valu = valu + (gradu(idim,:,:) .* msh_side.normal(idim,:,:));
  end

  w = msh_side.quad_weights .* msh_side.jacdet;
  err_elem = sum (reshape (valu, [msh_side.nqn, msh_side.nel]) .* w);
  err = sum (err_elem);

  flux = flux + err;
end

end
