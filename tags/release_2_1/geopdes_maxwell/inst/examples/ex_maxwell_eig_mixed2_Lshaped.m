% EX_MAXWELL_EIG_MIXED2_LSHAPED: solve Maxwell eigenproblem in the L-shaped domain, with the second mixed formulation.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_Lshaped_C0.txt';
%problem_data.geo_name = 'geo_Lshaped_C1.txt';

% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];

% Physical parameters
problem_data.c_elec_perm = @(x, y) ones(size(x));
problem_data.c_magn_perm = @(x, y) ones(size(x));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data 
method_data.degree     = [3 3];     % Degree of the bsplines
method_data.regularity = [2 2];     % Regularity of the splines
method_data.nsub       = [8 8];     % Number of subdivisions
method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, sp_mul, eigv, eigf] = ...
                      solve_maxwell_eig_mixed2_2d (problem_data, method_data);

% 4) POSTPROCESSING
[eigv, perm] = sort (eigv);

fprintf ('First computed eigenvalues: \n')
disp (eigv(1:6))

vtk_pts = {linspace(0,1,30) linspace(0,1,30)};
figure
[eu, F] = sp_eval (eigf(:,perm(9)), space, geometry, vtk_pts);
quiver (squeeze(F(1,:,:)), squeeze(F(2,:,:)), ...
        squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
title ('8^{th} eigenfunction')

%!demo
%! ex_maxwell_eig_mixed2_Lshaped

%!test
%! problem_data.geo_name = 'geo_Lshaped_C0.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4];
%! problem_data.c_elec_perm = @(x, y) ones(size(x));
%! problem_data.c_magn_perm = @(x, y) ones(size(x));
%! method_data.degree     = [3 3];     % Degree of the bsplines
%! method_data.regularity = [2 2];     % Regularity of the splines
%! method_data.nsub       = [8 8];     % Number of subdivisions
%! method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule
%! [geometry, msh, space, sp_mul, eigv, eigf] = solve_maxwell_eig_mixed2_2d (problem_data, method_data);
%! [eigv, perm] = sort (eigv);
%! assert (msh.nel, 128)
%! assert (space.ndof, 430)
%! assert (sp_mul.ndof, 200)
%! assert (eigv(2:6), [1.47283540635722; 3.53399368976646; 9.86962050365007; 9.86962103171719; 11.38943951321592], 1e-13)
