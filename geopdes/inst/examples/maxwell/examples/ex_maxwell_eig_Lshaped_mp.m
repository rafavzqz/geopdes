% EX_MAXWELL_EIG_LSHAPED_MP: solve Maxwell eigenproblem in the multipatch L-shaped domain.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_Lshaped_mp.txt';
%problem_data.geo_name = 'geo_Lshaped_mp_b.txt';
%problem_data.geo_name = 'geo_Lshaped_mp_c.txt';

% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4 5 6];

% Physical parameters
problem_data.c_elec_perm = @(x, y) ones(size(x));
problem_data.c_magn_perm = @(x, y) ones(size(x));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data 
method_data.degree     = [3 3];     % Degree of the bsplines
method_data.regularity = [2 2];     % Regularity of the splines
method_data.nsub       = [6 6];     % Number of subdivisions
method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, eigv, eigf] = mp_solve_maxwell_eig (problem_data, method_data);

% 4) POSTPROCESSING
[eigv, perm] = sort (eigv);
nzeros = numel (find (eigv < 1e-10));

fprintf ('Number of zero eigenvalues: %i \n', nzeros)
fprintf ('First nonzero eigenvalues: \n')
disp (eigv(nzeros+1:nzeros+5))

%!demo
%! ex_maxwell_eig_Lshaped_mp

%!test
%! problem_data.geo_name = 'geo_Lshaped_mp.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_elec_perm = @(x, y) ones(size(x));
%! problem_data.c_magn_perm = @(x, y) ones(size(x));
%! method_data.degree     = [3 3];     % Degree of the bsplines
%! method_data.regularity = [2 2];     % Regularity of the splines
%! method_data.nsub       = [6 6];     % Number of subdivisions
%! method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule
%! [geometry, msh, space, eigv, eigf] = mp_solve_maxwell_eig (problem_data, method_data);
%! [eigv, perm] = sort (eigv);
%! nzeros = numel (find (eigv < 1e-10));
%! assert (msh.nel, 108)
%! assert (space.ndof, 416)
%! assert (nzeros, 161)
%! assert (eigv(nzeros+(1:5)), [1.47383596756687; 3.53401518183127; 9.86961195350075; 9.86961195350077; 11.38946326693307], 1e-13)

%!test
%! problem_data.geo_name = 'geo_Lshaped_mp_b.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_elec_perm = @(x, y) ones(size(x));
%! problem_data.c_magn_perm = @(x, y) ones(size(x));
%! method_data.degree     = [3 3];     % Degree of the bsplines
%! method_data.regularity = [2 2];     % Regularity of the splines
%! method_data.nsub       = [6 6];     % Number of subdivisions
%! method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule
%! [geometry, msh, space, eigv, eigf] = mp_solve_maxwell_eig (problem_data, method_data);
%! [eigv, perm] = sort (eigv);
%! nzeros = numel (find (eigv < 1e-10));
%! assert (msh.nel, 108)
%! assert (space.ndof, 416)
%! assert (nzeros, 161)
%! assert (eigv(nzeros+(1:5)), [1.47383596756687; 3.53401518183127; 9.86961195350075; 9.86961195350077; 11.38946326693307], 1e-13)

%!test
%! problem_data.geo_name = 'geo_Lshaped_mp_c.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_elec_perm = @(x, y) ones(size(x));
%! problem_data.c_magn_perm = @(x, y) ones(size(x));
%! method_data.degree     = [3 3];     % Degree of the bsplines
%! method_data.regularity = [2 2];     % Regularity of the splines
%! method_data.nsub       = [6 6];     % Number of subdivisions
%! method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule
%! [geometry, msh, space, eigv, eigf] = mp_solve_maxwell_eig (problem_data, method_data);
%! [eigv, perm] = sort (eigv);
%! nzeros = numel (find (eigv < 1e-10));
%! assert (msh.nel, 108)
%! assert (space.ndof, 416)
%! assert (nzeros, 161)
%! assert (eigv(nzeros+(1:5)), [1.47383596756687; 3.53401518183127; 9.86961195350075; 9.86961195350077; 11.38946326693307], 1e-13)