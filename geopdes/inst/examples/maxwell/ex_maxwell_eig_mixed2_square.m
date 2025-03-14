% EX_MAXWELL_EIG_MIXED2_SQUARE: solve Maxwell eigenproblem in the unit square, with the second mixed formulation.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_square.txt';

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
method_data.nsub       = [10 10];     % Number of subdivisions
method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, sp_mul, eigv, eigf] = ...
                       solve_maxwell_eig_mixed2_2d (problem_data, method_data);

% 4) POSTPROCESSING
[eigv, perm] = sort (abs(eigv));

fprintf ('First computed eigenvalues: \n')
disp (eigv(1:6))

figure
sp_plot_solution (eigf(1:space.ndof,perm(9)), space, geometry, [30 30])
title ('8^{th} eigenfunction')

%!demo
%! ex_maxwell_eig_mixed2_square

%!test
%! problem_data.geo_name = 'geo_square.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4];
%! problem_data.c_elec_perm = @(x, y) ones(size(x));
%! problem_data.c_magn_perm = @(x, y) ones(size(x));
%! method_data.degree     = [3 3];     % Degree of the bsplines
%! method_data.regularity = [2 2];     % Regularity of the splines
%! method_data.nsub       = [10 10];     % Number of subdivisions
%! method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule
%! [geometry, msh, space, sp_mul, eigv, eigf] = solve_maxwell_eig_mixed2_2d (problem_data, method_data);
%! [eigv, perm] = sort (eigv);
%! eigv = eigv(eigv>0 & eigv<Inf);
%! nzeros = numel (find (eigv < 1e-4));
%! assert (msh.nel, 100)
%! assert (space.ndof, 312)
%! assert (sp_mul.ndof, 144)
%! assert (eigv(nzeros+(1:5))/pi^2, [1.00000003326399; 1.00000003326400; 2.00000006652799; 4.00000968384168; 4.00000968384168], 2e-14)
