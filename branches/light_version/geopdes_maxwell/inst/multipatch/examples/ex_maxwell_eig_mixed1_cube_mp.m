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
method_data.nsub       = [3 3 3]; % Number of subdivisions
method_data.nquad      = [3 3 3]; % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, sp_mul, eigv, eigf, gnum, dofs_ornt, gnum_mul] = ...
                     mp_solve_maxwell_eig_mixed1 (problem_data, method_data);

% 4) POSTPROCESSING
[eigv, perm] = sort (eigv);

fprintf ('First computed eigenvalues: \n')
disp (eigv(1:5))

%!demo
%! ex_maxwell_eig_mixed1_cube_mp

%!test
%! problem_data.geo_name = 'geo_2cubesa.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_elec_perm = @(x, y, z) ones(size(x));
%! problem_data.c_magn_perm = @(x, y, z) ones(size(x));
%! method_data.degree     = [2 2 2]; % Degree of the bsplines
%! method_data.regularity = [1 1 1]; % Regularity of the splines
%! method_data.nsub       = [3 3 3]; % Number of subdivisions
%! method_data.nquad      = [3 3 3]; % Points for the Gaussian quadrature rule
%! [geometry, msh, space, sp_mul, eigv, eigf, gnum, dofs_ornt, gnum_mul] = mp_solve_maxwell_eig_mixed1 (problem_data, method_data);
%! [eigv, perm] = sort (eigv);
%! assert (sum (cellfun (@(x) x.nel, msh)), 108)
%! assert (max([gnum{:}]), 923)
%! assert (max([gnum_mul{:}]), 360)
%! assert (eigv(1:5), [19.7435520600124; 19.7607903571162; 19.7629814831161; 29.6336619501220; 29.6336619501221], 5e-12)

%!test
%! problem_data.geo_name = 'geo_2cubesb.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_elec_perm = @(x, y, z) ones(size(x));
%! problem_data.c_magn_perm = @(x, y, z) ones(size(x));
%! method_data.degree     = [2 2 2]; % Degree of the bsplines
%! method_data.regularity = [1 1 1]; % Regularity of the splines
%! method_data.nsub       = [3 3 3]; % Number of subdivisions
%! method_data.nquad      = [3 3 3]; % Points for the Gaussian quadrature rule
%! [geometry, msh, space, sp_mul, eigv, eigf, gnum, dofs_ornt, gnum_mul] = mp_solve_maxwell_eig_mixed1 (problem_data, method_data);
%! [eigv, perm] = sort (eigv);
%! assert (sum (cellfun (@(x) x.nel, msh)), 108)
%! assert (max([gnum{:}]), 923)
%! assert (max([gnum_mul{:}]), 360)
%! assert (eigv(1:5), [19.7435520600124; 19.7607903571162; 19.7629814831161; 29.6336619501220; 29.6336619501221], 5e-12)

%!test
%! problem_data.geo_name = 'geo_2cubesc.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_elec_perm = @(x, y, z) ones(size(x));
%! problem_data.c_magn_perm = @(x, y, z) ones(size(x));
%! method_data.degree     = [2 2 2]; % Degree of the bsplines
%! method_data.regularity = [1 1 1]; % Regularity of the splines
%! method_data.nsub       = [3 3 3]; % Number of subdivisions
%! method_data.nquad      = [3 3 3]; % Points for the Gaussian quadrature rule
%! [geometry, msh, space, sp_mul, eigv, eigf, gnum, dofs_ornt, gnum_mul] = mp_solve_maxwell_eig_mixed1 (problem_data, method_data);
%! [eigv, perm] = sort (eigv);
%! assert (sum (cellfun (@(x) x.nel, msh)), 108)
%! assert (max([gnum{:}]), 923)
%! assert (max([gnum_mul{:}]), 360)
%! assert (eigv(1:5), [19.7435520600124; 19.7607903571162; 19.7629814831161; 29.6336619501220; 29.6336619501221], 5e-12)

%!test
%! problem_data.geo_name = 'geo_2cubesd.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_elec_perm = @(x, y, z) ones(size(x));
%! problem_data.c_magn_perm = @(x, y, z) ones(size(x));
%! method_data.degree     = [2 2 2]; % Degree of the bsplines
%! method_data.regularity = [1 1 1]; % Regularity of the splines
%! method_data.nsub       = [3 3 3]; % Number of subdivisions
%! method_data.nquad      = [3 3 3]; % Points for the Gaussian quadrature rule
%! [geometry, msh, space, sp_mul, eigv, eigf, gnum, dofs_ornt, gnum_mul] = mp_solve_maxwell_eig_mixed1 (problem_data, method_data);
%! [eigv, perm] = sort (eigv);
%! assert (sum (cellfun (@(x) x.nel, msh)), 108)
%! assert (max([gnum{:}]), 923)
%! assert (max([gnum_mul{:}]), 360)
%! assert (eigv(1:5), [19.7435520600124; 19.7607903571162; 19.7629814831161; 29.6336619501220; 29.6336619501221], 5e-12)

%!test
%! problem_data.geo_name = 'geo_2cubese.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_elec_perm = @(x, y, z) ones(size(x));
%! problem_data.c_magn_perm = @(x, y, z) ones(size(x));
%! method_data.degree     = [2 2 2]; % Degree of the bsplines
%! method_data.regularity = [1 1 1]; % Regularity of the splines
%! method_data.nsub       = [3 3 3]; % Number of subdivisions
%! method_data.nquad      = [3 3 3]; % Points for the Gaussian quadrature rule
%! [geometry, msh, space, sp_mul, eigv, eigf, gnum, dofs_ornt, gnum_mul] = mp_solve_maxwell_eig_mixed1 (problem_data, method_data);
%! [eigv, perm] = sort (eigv);
%! assert (sum (cellfun (@(x) x.nel, msh)), 108)
%! assert (max([gnum{:}]), 923)
%! assert (max([gnum_mul{:}]), 360)
%! assert (eigv(1:5), [19.7435520600124; 19.7607903571162; 19.7629814831161; 29.6336619501220; 29.6336619501221], 5e-12)

%!test
%! problem_data.geo_name = 'geo_2cubesf.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_elec_perm = @(x, y, z) ones(size(x));
%! problem_data.c_magn_perm = @(x, y, z) ones(size(x));
%! method_data.degree     = [2 2 2]; % Degree of the bsplines
%! method_data.regularity = [1 1 1]; % Regularity of the splines
%! method_data.nsub       = [3 3 3]; % Number of subdivisions
%! method_data.nquad      = [3 3 3]; % Points for the Gaussian quadrature rule
%! [geometry, msh, space, sp_mul, eigv, eigf, gnum, dofs_ornt, gnum_mul] = mp_solve_maxwell_eig_mixed1 (problem_data, method_data);
%! [eigv, perm] = sort (eigv);
%! assert (sum (cellfun (@(x) x.nel, msh)), 108)
%! assert (max([gnum{:}]), 923)
%! assert (max([gnum_mul{:}]), 360)
%! assert (eigv(1:5), [19.7435520600124; 19.7607903571162; 19.7629814831161; 29.6336619501220; 29.6336619501221], 5e-12)

%!test
%! problem_data.geo_name = 'geo_2cubesg.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_elec_perm = @(x, y, z) ones(size(x));
%! problem_data.c_magn_perm = @(x, y, z) ones(size(x));
%! method_data.degree     = [2 2 2]; % Degree of the bsplines
%! method_data.regularity = [1 1 1]; % Regularity of the splines
%! method_data.nsub       = [3 3 3]; % Number of subdivisions
%! method_data.nquad      = [3 3 3]; % Points for the Gaussian quadrature rule
%! [geometry, msh, space, sp_mul, eigv, eigf, gnum, dofs_ornt, gnum_mul] = mp_solve_maxwell_eig_mixed1 (problem_data, method_data);
%! [eigv, perm] = sort (eigv);
%! assert (sum (cellfun (@(x) x.nel, msh)), 108)
%! assert (max([gnum{:}]), 923)
%! assert (max([gnum_mul{:}]), 360)
%! assert (eigv(1:5), [19.7435520600124; 19.7607903571162; 19.7629814831161; 29.6336619501220; 29.6336619501221], 5e-12)

%!test
%! problem_data.geo_name = 'geo_2cubesh.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_elec_perm = @(x, y, z) ones(size(x));
%! problem_data.c_magn_perm = @(x, y, z) ones(size(x));
%! method_data.degree     = [2 2 2]; % Degree of the bsplines
%! method_data.regularity = [1 1 1]; % Regularity of the splines
%! method_data.nsub       = [3 3 3]; % Number of subdivisions
%! method_data.nquad      = [3 3 3]; % Points for the Gaussian quadrature rule
%! [geometry, msh, space, sp_mul, eigv, eigf, gnum, dofs_ornt, gnum_mul] = mp_solve_maxwell_eig_mixed1 (problem_data, method_data);
%! [eigv, perm] = sort (eigv);
%! assert (sum (cellfun (@(x) x.nel, msh)), 108)
%! assert (max([gnum{:}]), 923)
%! assert (max([gnum_mul{:}]), 360)
%! assert (eigv(1:5), [19.7435520600124; 19.7607903571162; 19.7629814831161; 29.6336619501220; 29.6336619501221], 5e-12)
