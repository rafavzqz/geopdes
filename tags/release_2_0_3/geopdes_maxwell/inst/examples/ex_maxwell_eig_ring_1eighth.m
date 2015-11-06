% EX_MAXWELL_EIG_RING_1EIGHTH: solve Maxwell eigenproblem in 1/8 of a ring.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_ring_1eighth.txt';

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
[geometry, msh, space, eigv, eigf] = ...
                             solve_maxwell_eig_2d (problem_data, method_data);

% 4) POSTPROCESSING
[eigv, perm] = sort (eigv);
nzeros = numel (find (eigv < 1e-10));

fprintf ('Number of zero eigenvalues: %i \n', nzeros)
fprintf ('First nonzero eigenvalues: \n')
disp (eigv(nzeros+1:nzeros+5))

vtk_pts = {linspace(0,1,30) linspace(0,1,30)};
figure
[eu, F] = sp_eval (eigf(:,perm(nzeros+8)), space, geometry, vtk_pts);
quiver (squeeze(F(1,:,:)), squeeze(F(2,:,:)), ...
        squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
title ('8^{th} eigenfunction')

%!demo
%! ex_maxwell_eig_ring_1eighth
