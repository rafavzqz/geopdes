% SOLVE_CAHN_HILLIARD: solve the Cahn-Hilliard equation, with a generalized alpha discretization in time.
%
% The functions solves the problem of finding u such that
%
%  du/dt - Delta (mu(u) - lambda*Delta u) = 0
%
% with Delta the Laplacian, and mu(u) = alpha u^3 - beta u, and periodic boundary conditions.
%
% The values of alpha and beta (or mu itself) can be changed in op_gradmu_gradv_tp.
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
%    - Time_max:     final time
%    - fun_u:        initial condition. Equal to zero by default.
%    - fun_udot:     initial condition for time derivative. Equal to zero by default.
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%    - dt:         time step size for generalized-alpha method
%    - rho_inf_gen_alpha: parameter in [0,1], which governs numerical damping of the generalized alpha method
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
% Copyright (C) 2023 Michele Torre, Rafael Vazquez
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

function [geometry, msh, space, results] = solve_cahn_hilliard_mpC1 (problem_data, method_data, save_info)

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




% Construct geometry structure, and information for interfaces and boundaries
[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (geo_name);
npatch = numel (geometry);

msh = cell (1, npatch); 
sp = cell (1, npatch);
for iptc = 1:npatch

% Define the refined mesh, with tensor product structure
  [knots{iptc}, zeta{iptc}] = ...
         kntrefine (geometry(iptc).nurbs.knots, nsub-1, degree, regularity);

% Compute the quadrature rule
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (zeta{iptc}, rule);
  msh{iptc} = msh_cartesian (zeta{iptc}, qn, qw, geometry(iptc));

% Evaluate the discrete space basis functions in the quadrature points
  sp{iptc} = sp_bspline (knots{iptc}, degree, msh{iptc});
end

[edges, vertices] = vertices_struct (geometry, interfaces, boundaries, boundary_interfaces);
msh = msh_multipatch (msh, boundaries);
% space = sp_multipatch (sp, msh, interfaces, boundary_interfaces);
space = sp_multipatch_C1 (sp, msh, geometry, edges, vertices);
clear sp

% % Compute and assemble the matrices 
% c_diff  = @(x, y)  ones(size(x));
% stiff_mat = op_gradu_gradv_mp (space, space, msh, c_diff);
% f  = @(x, y)  ones(size(x));
% rhs = op_f_v_mp (space, msh, f);
% drchlt_dofs = 1:10;
% int_dofs = setdiff (1:space.ndof, drchlt_dofs);
% % Solve the linear system
% u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);
% output_file = 'Lshaped_mp_BSP_Deg3_Reg2_Sub9';
% vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
% fprintf ('The result is saved in the file %s.pvd \n \n', output_file);
% sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')
% figure
% sp_plot_solution (u, space, geometry, vtk_pts)
% return


%%-------------------------------------------------------------------------
% Generalized alpha parameters
a_m = .5*((3-rho_inf_gen_alpha)/(1+rho_inf_gen_alpha));
a_f =  1/(1+rho_inf_gen_alpha);
gamma =  .5 + a_m - a_f;

%%-------------------------------------------------------------------------
% No flux b.c. (essential boundary condition)
% Set Neumann boundary conditions for non-periodic sides
nmnn_bou   = 1:numel(boundaries);

%%-------------------------------------------------------------------------
% Precompute matrices

% Compute the mass matrix
mass_mat = op_u_v_mp(space,space,msh);

% Compute the laplace matrix
lapl_mat = op_laplaceu_laplacev_mp (space, space, msh, lambda);

% Compute the boundary term
term4 = int_term_4 (space, msh, lambda, nmnn_bou);

% Compute the penalty matrix
[Pen, pen_rhs] = penalty_matrix (space, msh, nmnn_bou, pen_nitsche);


%%-------------------------------------------------------------------------
% Initial conditions
time = 0;
if (exist('fun_u', 'var') && ~isempty(fun_u))
  rhs = op_f_v_mp (space, msh, fun_u);
  u_n = (mass_mat + pen_projection/pen_nitsche * Pen)\rhs;
else
  u_n = zeros(space.ndof, 1);
end

if (exist('fun_udot', 'var') && ~isempty(fun_udot))
  rhs = op_f_v_mp(space, msh, fun_udot);
  udot_n = (mass_mat + pen_projection/pen_nitsche * Pen)\rhs;
else
  udot_n = zeros(space.ndof, 1);
end

% flux_initial = check_flux_phase_field(space, msh, u_n, zeros(size(u_n)));
% disp(strcat('initial flux =',num2str(flux_initial)))

%%-------------------------------------------------------------------------
% Initialize structure to store the results
save_id = 1;
results.u = zeros(length(u_n), length(save_info)+1);
results.udot = zeros(length(u_n), length(save_info)+1);
results.time = zeros(length(save_info)+1,1);
flag_stop_save = false;

% Save initial conditions
results.u(:,1) = u_n;
results.udot(:,1) = udot_n;
results.time(1) = time;

%%-------------------------------------------------------------------------
% Loop over time steps
while time < Time_max
  disp('----------------------------------------------------------------')
  disp(strcat('time step t=',num2str(time)))

  [u_n1, udot_n1] = generalized_alpha_step(u_n, udot_n, dt, a_m, a_f, gamma, ...
                    mass_mat, lapl_mat, term4, Pen, pen_rhs, space, msh);

  % check flux through the boundary
%   flux = check_flux_phase_field(space, msh, u_n1, zeros(size(u_n1)));
%   disp(strcat('flux norm =',num2str(flux)))

  % Time step update
  time = time + dt;
  u_n = u_n1;
  udot_n = udot_n1;

  % Check max time
  if (time + dt > Time_max)
    dt = Time_max-time;
  end

  % Store results
  if (flag_stop_save == false)
    if (time >= save_info(save_id))
      save_id = save_id + 1;
      results.u(:,save_id) = u_n;
      results.udot(:,save_id) = udot_n;
      results.time(save_id) = time;
      if (save_id > length(save_info))
        flag_stop_save = true;
      end
    end
  end
end

disp('----------------------------------------------------------------')
disp(strcat('END ANALYSIS t=',num2str(time)))
disp('----------------------------------------------------------------')

% Crop results
results.u = results.u(:,1:save_id);
results.udot = results.udot(:,1:save_id);
results.time = results.time(1:save_id);

end

%--------------------------------------------------------------------------
% FUNCTIONS
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% One step of generalized alpha-method
%--------------------------------------------------------------------------

function [u_n1, udot_n1] = generalized_alpha_step(u_n, udot_n, dt, a_m, a_f, gamma, ...
                       mass_mat, lapl_mat, term4, Pen, pen_rhs, space, msh)

% Convergence criteria
  n_max_iter = 20;
  tol_rel_res = 1e-10;
  tol_abs_res = 1e-10;

% Predictor step
  u_n1 = u_n;
  udot_n1 = (gamma-1)/gamma * udot_n; 

% Newton loop
  for iter = 0:n_max_iter

  % Field at alpha level
    udot_a = udot_n + a_m *(udot_n1-udot_n);
    u_a = u_n + a_f *(u_n1-u_n);

  % Compute the residual (internal)
    [Res_gl, stiff_mat] = Res_K_cahn_hilliard (space, msh, ...
                          mass_mat, lapl_mat, term4, Pen, pen_rhs, ...
                          u_a, udot_a);

  % Convergence check
    if iter == 0
      norm_res_0 = norm(Res_gl);
    end
    norm_res = norm(Res_gl);
    

    if norm_res/norm_res_0 < tol_rel_res
      disp(strcat('iteration n°=',num2str(iter)))
      disp(strcat('norm (abs) residual=',num2str(norm_res)))
      break
    end
    if norm_res<tol_abs_res
      disp(strcat('iteration n°=',num2str(iter)))
      disp(strcat('norm absolute residual=',num2str(norm_res)))
      break
    end
    if iter == n_max_iter
      disp(strcat('Newton reached the maximum number of iterations'))
      disp(strcat('norm residual=',num2str(norm_res)))
    end
    
  % Compute the update, and update the solution
    A_gl = a_m * mass_mat + a_f * gamma *dt * stiff_mat ; 
    d_udot = - A_gl\Res_gl;

    udot_n1 = udot_n1 + d_udot;
    u_n1 = u_n1 + gamma * dt* d_udot;
end

end

%--------------------------------------------------------------------------
% Canh-Hilliard residual and tangent matrix
%--------------------------------------------------------------------------

function [Res_gl, stiff_mat] = Res_K_cahn_hilliard(space, msh, ...
                          mass_mat, lapl_mat, term4, Pen, pen_rhs, ...
                          u_a, udot_a)

    % Double well (matrices)
    [term2, term2K] = op_gradmu_gradv_mp(space, msh, u_a);    
 
    % Residual
    Res_gl = mass_mat*udot_a + term2*u_a  + lapl_mat*u_a;

    % Tangent stiffness matrix (mass is not considered here)
    stiff_mat = term2 + term2K + lapl_mat;

    % in case of neumann BC, add boundary terms
    if isempty(term4) == 0
         Res_gl = Res_gl - (term4 + term4') * u_a + Pen*u_a - pen_rhs;
         stiff_mat = stiff_mat - (term4 + term4') + Pen;
    end
end

%--------------------------------------------------------------------------
% Integral of the double-well function
%--------------------------------------------------------------------------

function [A, B] = op_gradmu_gradv_mp (space, msh,  uhat)



  A = spalloc (space.ndof, space.ndof, 3*space.ndof);
  B = spalloc (space.ndof, space.ndof, 6*space.ndof);



  for iptc = 1:msh.npatch
    [Cpatch_u, Cpatch_cols_u] = sp_compute_Cpatch (space, iptc);
    u_ptc = Cpatch_u * uhat(Cpatch_cols_u);


    [Ap, Bp] = op_gradmu_gradv_tp (space.sp_patch{iptc}, msh.msh_patch{iptc}, u_ptc );

    
    
    A(Cpatch_cols_u,Cpatch_cols_u) = ...
          A(Cpatch_cols_u,Cpatch_cols_u) + Cpatch_u.' * Ap * Cpatch_u;

    B(Cpatch_cols_u,Cpatch_cols_u) = ...
          B(Cpatch_cols_u,Cpatch_cols_u) + Cpatch_u.' * Bp * Cpatch_u;
  end


end


function [A, B] = op_gradmu_gradv_tp (space, msh,  uhat)

% Coefficients of the double well function.
  alpha = 1;
  beta = 1;
  A = spalloc (space.ndof, space.ndof, 3*space.ndof);
  B = spalloc (space.ndof, space.ndof, 6*space.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col, 'gradient', true);

    % Evaluate the field and its gradient at the Gaussian points
    utemp = sp_eval_msh (uhat, sp_col, msh_col, {'value', 'gradient'});
    u = utemp{1};
    gradu = utemp{2};

    % Polynomial formulation for the double-well
    coeffs_A = 3.* alpha .* u.^2 - beta;
    A = A + op_gradu_gradv (sp_col, sp_col, msh_col, coeffs_A);

    coeffs_B = 6.* alpha .* u;
    coeffs_Bv = gradu;
    for idim = 1:msh.ndim
      coeffs_Bv(idim,:,:) = coeffs_Bv(idim,:,:) .* reshape(coeffs_B, 1, size(coeffs_B,1), size(coeffs_B,2));
    end
    B = B + op_vel_dot_gradu_v (sp_col, sp_col, msh_col, coeffs_Bv)';
  end

end



%--------------------------------------------------------------------------
% term 4, boundary term
%--------------------------------------------------------------------------

function [A] = int_term_4 (space, msh,  lambda, nmnn_sides)

  if length(nmnn_sides)>0

    A =  spalloc (space.ndof, space.ndof, 3*space.ndof);

    for iref = nmnn_sides
    for bnd_side = 1:msh.boundaries(iref).nsides  
        iptc = msh.boundaries(iref).patches(bnd_side);
        iside = msh.boundaries(iref).faces(bnd_side);

        msh_side = msh_eval_boundary_side (msh.msh_patch{iptc}, iside ) ;
        msh_side_int = msh_boundary_side_from_interior (msh.msh_patch{iptc}, iside ) ;

        sp_side = space.sp_patch{iptc}.constructor (msh_side_int);
        sp_side = sp_precompute ( sp_side , msh_side_int , 'gradient' , true, 'laplacian', true );

        for idim = 1:msh.rdim
            x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
        end
        coe_side = lambda (x{:});
        tmp = op_gradv_n_laplaceu(sp_side ,sp_side ,msh_side, coe_side);


        [Cpatch, Cpatch_cols] = sp_compute_Cpatch (space, iptc);
        A(Cpatch_cols,Cpatch_cols) = ...
              A(Cpatch_cols,Cpatch_cols) + Cpatch.' * tmp * Cpatch;

    end
    end

else
    A = [];
end


end

%--------------------------------------------------------------------------
% Operator grav_n_laplaceu
%--------------------------------------------------------------------------

function varargout = op_gradv_n_laplaceu (spu, spv, msh, coeff)

  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], ...
		   msh.nqn, spv.nsh_max, msh.nel);

  ndim = size (gradv, 2);

  shpu = reshape (spu.shape_function_laplacians, spu.ncomp, msh.nqn, spu.nsh_max, msh.nel);

  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  jacdet_weights = msh.jacdet .* msh.quad_weights .* coeff;
  
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
      gradv_iel = gradv(:,:,:,:,iel);
      normal_iel = reshape (msh.normal(:,:,iel), [1, ndim, msh.nqn]);

      gradv_n = reshape (sum (bsxfun (@times, gradv_iel, normal_iel), 2), spv.ncomp, msh.nqn, spv.nsh_max, 1);
      shpu_iel = reshape (shpu(:, :, :, iel), spu.ncomp, msh.nqn, 1, spu.nsh_max);

      jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);

      gradv_n_times_jw = bsxfun (@times, jacdet_iel, gradv_n);
      tmp1 = sum (bsxfun (@times, gradv_n_times_jw, shpu_iel), 1);
      elementary_values = reshape (sum (tmp1, 2), spv.nsh_max, spu.nsh_max);
      
      [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spu.connectivity(:,iel));
      indices = rows_loc & cols_loc;
      rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = rows_loc(indices);
      cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = cols_loc(indices);
      values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = elementary_values(indices);
      ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel);
      
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradv_n_u: singular map in element number %d', iel)
    end
  end
  

  if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows, cols, values, spv.ndof, spu.ndof);
  elseif (nargout == 3)
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = values;
  else
    error ('op_gradv_n_u: wrong number of output arguments')
  end
  
  
end

%--------------------------------------------------------------------------
% penalty term
%--------------------------------------------------------------------------

function [P, rhs] = penalty_matrix (space, msh, nmnn_sides, pen)

    P =  spalloc (space.ndof, space.ndof, 3*space.ndof);
    rhs = zeros(space.ndof,1);

    for iref = nmnn_sides
    for bnd_side = 1:msh.boundaries(iref).nsides
      iptc = msh.boundaries(iref).patches(bnd_side);
      iside = msh.boundaries(iref).faces(bnd_side);

      msh_side = msh_eval_boundary_side (msh.msh_patch{iptc}, iside);
      msh_side_from_interior = msh_boundary_side_from_interior (msh.msh_patch{iptc}, iside);

      sp_bnd = space.sp_patch{iptc}.constructor (msh_side_from_interior);
      sp_bnd = sp_precompute (sp_bnd, msh_side_from_interior, 'value', false, 'gradient', true);

      coe_side = pen .* msh_side.charlen; 
      tmp = op_gradu_n_gradv_n(sp_bnd, sp_bnd, msh_side, coe_side);
      
      [Cpatch, Cpatch_cols] = sp_compute_Cpatch (space, iptc);
      P(Cpatch_cols,Cpatch_cols) = P(Cpatch_cols,Cpatch_cols) + ...
        Cpatch.' * (tmp) * Cpatch;


      %rhs = rhs + rhs_pen;

    end
    end

end


% function [mass_pen, rhs_pen] = penalty_grad (space, msh, i_side, pen)
% 
% msh_side = msh_eval_boundary_side (msh, i_side ) ;
% msh_side_int = msh_boundary_side_from_interior (msh, i_side ) ;
% sp_side = space.constructor ( msh_side_int) ;
% sp_side = sp_precompute ( sp_side , msh_side_int , 'gradient', true );
% 
% 
% coe_side = pen .* msh_side.charlen; 
% mass_pen = op_gradu_n_gradv_n(sp_side, sp_side, msh_side, coe_side);
% 
% 
% rhs_pen  = zeros(space.ndof,1); % no flux
% 
% end

%--------------------------------------------------------------------------
% check flux through the boundaries
%--------------------------------------------------------------------------

% function flux = check_flux_phase_field(space, msh, uhat, uhat0)
% 
% sides = [1,2,3,4]; % all the boundaries
% flux = 0;
% 
% for iside=1:length(sides)   
% 
%     msh_side = msh_eval_boundary_side (msh, sides(iside) ) ;
%     msh_side_int = msh_boundary_side_from_interior (msh, sides(iside) ) ;
%     sp_side = space.constructor ( msh_side_int) ;
%     sp_side = sp_precompute ( sp_side , msh_side_int , 'gradient' , true );
% 
% 
%     gradu = sp_eval_msh (uhat-uhat0, sp_side, msh_side, 'gradient');
%    
% 
%     
%     valu = zeros(sp_side.ncomp, size(msh_side.quad_weights,1), size(msh_side.quad_weights,2));
%     for idim = 1:msh.rdim
%         valu = valu + (gradu(idim,:,:) .* msh_side.normal(idim,:,:));
%     end
% 
% 
%     w =msh_side.quad_weights .* msh_side.jacdet;
%     err_elem = sum (reshape (valu, [msh_side.nqn, msh_side.nel]) .* w);
%     err  = sum (err_elem);
% 
% 
%     flux = flux + err;
% 
% end
% 
% 
% 
% end

