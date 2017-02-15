% EX_COLLOCATION_LAPLACE_RING: solve the Laplace problem with a NURBS discretization by isogeometric collocation

% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_ring.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];
        
% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));
        
% Source and boundary terms
problem_data.f = @(x,y) 2*x.*(22.*x.^2.*y.^2+21.*y.^4-45.*y.^2+x.^4-5.*x.^2+4);
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) zeros(size(x));
        
% Exact solution (optional)
problem_data.uex     = @(x, y) -(x.^2+y.^2-1).*(x.^2+y.^2.-4).*x.*y.^2;
problem_data.graduex = @(x, y) cat (1,  reshape (-2*(x.*y).^2.*((x.^2+y.^2-1)+(x.^2+y.^2-4)) - (x.^2+y.^2-1).*(x.^2+y.^2-4).*y.^2, [1, size(x)]), ...
  reshape ( -2*x.*y.^3.*((x.^2+y.^2-1)+(x.^2+y.^2-4)) -   2*x.*y.*(x.^2+y.^2-1).*(x.^2+y.^2-4), [1, size(x)]));
        
% Discretization parameters (p and h)
method_data.degree     = [4 4]; % Degree of the splines in each direction
method_data.regularity = [3 3]; % Regularity of the splines, should be at least C^1.
method_data.nsub       = [8 8]; % Divide each subinterval of the original knot span in nsub subintervals
method_data.pts_case   = 2;     % Collocation points. 1: uniform, 2: Greville abscissae
% method_data.ncoll_pts  = [40 40]; % Number of collocation points in each direction. Only used for pts_case = 1.

% Solve with collocation method
[geometry, msh_coll, space_coll, u_coll] = solve_laplace_collocation (problem_data, method_data);

% Solve with Galerkin, for comparison
method_data.nquad = method_data.degree + 1;
[~, msh_gal, space_gal, u_gal] = solve_laplace_iso (problem_data, method_data);

% Plot of solution
vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
figure; 
[eu, F] = sp_eval (u_coll, space_coll, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
subplot (1,3,1)
surf (X, Y, eu)
title ('Collocation solution'), axis tight
[eu, F] = sp_eval (u_coll, space_coll, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
subplot (1,3,2)
surf (X, Y, eu)
title ('Galerkin solution'), axis tight
subplot (1,3,3)
surf (X, Y, problem_data.uex (X,Y))
title ('Exact solution'), axis tight

% Compute errors of collocation and Galerkin. Since it involves quadrature,
%  it is computed using the mesh and space objects from the Galerkin method.
disp ('Error in H1 and L2 norms, for isogeometric collocation')
[error_h1_coll, error_l2_coll] = sp_h1_error (space_gal, msh_gal, u_coll, problem_data.uex, problem_data.graduex);
disp([error_l2_coll error_h1_coll])

disp ('Error in H1 and L2 norms, for isogeometric Galerkin')
[error_h1_gal, error_l2_gal] = sp_h1_error (space_gal, msh_gal, u_gal, problem_data.uex, problem_data.graduex);
disp ([error_l2_gal error_h1_gal])
