% EX_MAXWELL_EIG_THICK_L: solve Maxwell eigenproblem in the thick L-shaped domain, with the first mixed formulation.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_thickL_C0.txt';
%problem_data.geo_name = 'geo_thickL_C1.txt';

% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4 5 6];

% Physical parameters
problem_data.c_elec_perm = @(x, y, z) ones(size(x));
problem_data.c_magn_perm = @(x, y, z) ones(size(x));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data 
method_data.degree     = [2 2 2]; % Degree of the bsplines
method_data.regularity = [1 1 1]; % Regularity of the splines
method_data.nsub       = [3 3 3]; % Number of subdivisions
method_data.nquad      = [3 3 3]; % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, sp_mul, eigv, eigf] = ...
                       solve_maxwell_eig_mixed1 (problem_data, method_data);

% 4) POSTPROCESSING
[eigv, perm] = sort (eigv);

fprintf ('First computed eigenvalues: \n')
disp (eigv(1:5))

%!demo
%! ex_maxwell_eig_mixed1_thick_L

%!test
%! problem_data.geo_name = 'geo_thickL_C0.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_elec_perm = @(x, y, z) ones(size(x));
%! problem_data.c_magn_perm = @(x, y, z) ones(size(x));
%! method_data.degree     = [2 2 2]; % Degree of the bsplines
%! method_data.regularity = [1 1 1]; % Regularity of the splines
%! method_data.nsub       = [3 3 3]; % Number of subdivisions
%! method_data.nquad      = [3 3 3]; % Points for the Gaussian quadrature rule
%! [geometry, msh, space, sp_mul, eigv, eigf] = solve_maxwell_eig_mixed1 (problem_data, method_data);
%! [eigv, perm] = sort (eigv);
%! assert (msh.nel, 54)
%! assert (space.ndof, 560)
%! assert (sp_mul.ndof, 225)
%! assert (eigv(1:5), [9.82837494805257; 11.34792479399617; 13.42675244735295; 15.29612032296900; 19.71848483816251], 1e-12)