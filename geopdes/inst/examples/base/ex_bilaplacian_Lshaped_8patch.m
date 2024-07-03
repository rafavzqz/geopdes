% PHYSICAL DATA OF THE PROBLEM
clear problem_data
problem_data.geo_name = 'geo_Lshaped_8patches.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4 5 6];
problem_data.weak_drchlt_sides = [];
        
% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));
% problem_data.grad_c_diff = @(x, y) cat (1, ...  
%             reshape (zeros(size(x)), [1, size(x)]), ...
%             reshape (zeros(size(x)), [1, size(x)]));

% Source and boundary terms
problem_data.f = @bilaplacian_rhs_Lshaped;
% problem_data.g = @(x, y, ind) bilaplacian_Lshaped_g_nmnn_8patches (x,y,ind);
problem_data.h = @(x, y, ind) solution_bilaplacian_Lshaped (x, y);

% Exact solution (optional)
problem_data.uex     = @(x, y) solution_bilaplacian_Lshaped (x, y);
problem_data.graduex = @(x, y) solution_bilaplacian_Lshaped_grad (x, y);
problem_data.hessuex = @(x, y) solution_bilaplacian_Lshaped_hessian (x, y);

% DISCRETIZATION PARAMETERS
clear method_data
method_data.degree      = [4 4];        % Degree of the splines
method_data.regularity  = [2 2];        % Regularity of the splines
method_data.nsub        = [16 16];      % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nquad       = method_data.degree+1; % Points for the Gaussian quadrature rule

[geometry, msh, space, u] = mp_solve_bilaplace_C1 (problem_data, method_data);

[err_h2,~,~,err_h2s] = sp_h2_error (space, msh, u, problem_data.uex, problem_data.graduex, problem_data.hessuex)

subplot(1,2,1)
npts = [51 51];
sp_plot_solution (u, space, geometry, npts); shading interp; title('Computed solution')
vtk_pts = {linspace(0,1,npts(1)), linspace(0,1,npts(2))};
subplot(1,2,2)
for iptc = 1:msh.npatch
  F = reshape (geometry(iptc).map(vtk_pts), [2, npts]);
  X = reshape (F(1,:), npts); Y = reshape (F(2,:), npts);
  surf(X, Y, problem_data.uex(X,Y));
  hold on; shading interp
end
title('Exact solution')
