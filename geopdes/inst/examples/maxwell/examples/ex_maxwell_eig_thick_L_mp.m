% EX_MAXWELL_EIG_THICK_L_MP: solve Maxwell eigenproblem in the multipatch thick L-shaped domain.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_thickL_mp.txt';
%problem_data.geo_name = 'geo_thickL_mp_b.txt';

% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4 5 6 7 8];

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
[geometry, msh, space, eigv, eigf] = mp_solve_maxwell_eig (problem_data, method_data);

% 4) POSTPROCESSING
[eigv, perm] = sort (eigv);
nzeros = numel (find (eigv < 1e-10));

fprintf ('Number of zero eigenvalues: %i \n', nzeros)
fprintf ('First nonzero eigenvalues: \n')
disp (eigv(nzeros+1:nzeros+5))

%!demo
%! ex_maxwell_eig_thick_L_mp

%!test
%! problem_data.geo_name = 'geo_thickL_mp.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6 7 8];
%! problem_data.c_elec_perm = @(x, y, z) ones(size(x));
%! problem_data.c_magn_perm = @(x, y, z) ones(size(x));
%! method_data.degree     = [2 2 2]; % Degree of the bsplines
%! method_data.regularity = [1 1 1]; % Regularity of the splines
%! method_data.nsub       = [3 3 3]; % Number of subdivisions
%! method_data.nquad      = [3 3 3]; % Points for the Gaussian quadrature rule
%! [geometry, msh, space, eigv, eigf] = mp_solve_maxwell_eig (problem_data, method_data);
%! [eigv, perm] = sort (eigv);
%! nzeros = numel (find (eigv < 1e-10));
%! assert (msh.nel, 81)
%! assert (space.ndof, 820)
%! assert (nzeros, 99)
%! assert (eigv(nzeros+(1:5)), [9.70264116281338; 11.35699119986159; 13.42422084425671; 15.24006039881071; 19.59275105292330], 1e-12)

%!test
%! problem_data.geo_name = 'geo_thickL_mp_b.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6 7 8];
%! problem_data.c_elec_perm = @(x, y, z) ones(size(x));
%! problem_data.c_magn_perm = @(x, y, z) ones(size(x));
%! method_data.degree     = [2 2 2]; % Degree of the bsplines
%! method_data.regularity = [1 1 1]; % Regularity of the splines
%! method_data.nsub       = [3 3 3]; % Number of subdivisions
%! method_data.nquad      = [3 3 3]; % Points for the Gaussian quadrature rule
%! [geometry, msh, space, eigv, eigf] = mp_solve_maxwell_eig (problem_data, method_data);
%! [eigv, perm] = sort (eigv);
%! nzeros = numel (find (eigv < 1e-10));
%! assert (msh.nel, 81)
%! assert (space.ndof, 820)
%! assert (nzeros, 99)
%! assert (eigv(nzeros+(1:5)), [9.70264116281338; 11.35699119986159; 13.42422084425671; 15.24006039881071; 19.59275105292330], 1e-12)