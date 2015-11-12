% EX_MAXWELL_EIG_MIXED1_FICHERA_MP: solve Maxwell eigenproblem in Fichera's corner domain, with the first mixed formulation.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_fichera.txt';

% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1];

% Physical parameters
problem_data.c_elec_perm = @(x, y, z) ones(size(x));
problem_data.c_magn_perm = @(x, y, z) ones(size(x));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data 
method_data.degree     = [2 2 2]; % Degree of the bsplines
method_data.regularity = [1 1 1]; % Regularity of the splines
method_data.nsub       = [2 2 2]; % Number of subdivisions
method_data.nquad      = [3 3 3]; % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, sp_mul, eigv, eigf] = ...
                  mp_solve_maxwell_eig_mixed1 (problem_data, method_data);

% 4) POSTPROCESSING
[eigv, perm] = sort (eigv);

fprintf ('First computed eigenvalues: \n')
disp (eigv(1:5))

%!demo
%! ex_maxwell_eig_mixed1_fichera_mp

%!test
%! problem_data.geo_name = 'geo_fichera.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1];
%! problem_data.c_elec_perm = @(x, y, z) ones(size(x));
%! problem_data.c_magn_perm = @(x, y, z) ones(size(x));
%! method_data.degree     = [2 2 2]; % Degree of the bsplines
%! method_data.regularity = [1 1 1]; % Regularity of the splines
%! method_data.nsub       = [2 2 2]; % Number of subdivisions
%! method_data.nquad      = [3 3 3]; % Points for the Gaussian quadrature rule
%! [geometry, msh, space, sp_mul, eigv, eigf] = mp_solve_maxwell_eig_mixed1 (problem_data, method_data);
%! [eigv, perm] = sort (eigv);
%! assert (msh.nel, 56)
%! assert (space.ndof, 801)
%! assert (sp_mul.ndof, 316)
%! assert (eigv(1:5), [3.16706451934816; 5.88666869517490; 5.88666869517490; 10.80766615203873; 10.82625504260329], 1e-12)
