% EX_MAXWELL_EIG_MIXED1_CUBE_MP: solve Maxwell eigenproblem in the unit cube with a multipatch domain and the first mixed formulation.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_2cubesa.txt';
% You can see how the patches should match in the files 
%  geo_2cubesa{b,c,d,e,f,g,h}.txt
% In all the eight cases the result must be the same

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
method_data.nsub       = [2 2 2]; % Number of subdivisions
method_data.nquad      = [3 3 3]; % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, sp_mul, eigv, eigf] = ...
                     mp_solve_maxwell_eig_mixed1_3d (problem_data, method_data);

% 4) POSTPROCESSING
[eigv, perm] = sort (eigv);

fprintf ('First computed eigenvalues: \n')
disp (eigv(1:5))

%!demo
%! ex_maxwell_eig_mixed1_cube_mp
