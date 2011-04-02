% EX_STOKES_2D_RT_SQUARE: solve the Stokes problem in the square using RT elements.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_square.txt';

% Type of boundary conditions for each side of the domain
problem_data.drchlt_sides = 1:4;
problem_data.nmnn_sides = [];

% Physical parameters
problem_data.viscosity = @(x, y) ones (size (x));

% Force term
fx = @(x, y) (exp(x).*(2.*(-1 + y).*y.*(2 + (-5 + y).*y) - ...
              8.*x.*(-1 + y).*y.*(2 + (-5 + y).*y) + ...
              6.*x.^3.*(-4 + (-5 + y).*(-2 + y).*y.*(1 + y)) + ...
              x.^2.*(12 + y.*(-38 + y.*(19 + (-6 + y).*y))) + ...
              x.^4.*(12 + y.*(-38 + y.*(19 + (-6 + y).*y)))));
fy = @(x, y) (456 - 912.*y + ...
              exp(x).*(-456 - 2.*x.*(-23 + 5.*x).*(10 + (-3 + x).*x) + ...
              2.*(456 + x.*(-466 + x.*(253 + x.*(-82 + 7.*x)))).* ...
              y  + (-6 + x.*(6 + x.*(-11 + x.*(22 + 7.*x)))).* ...
              y.^2 + 2.*(6 + x.*(10 + x.*(-29 + (-6 + x).*x))).* ...
              y.^3 + (-6 + x.*(3 + x).*(-2 + x.*(7 + x))).*y.^4));
problem_data.f  = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));

% Boundary terms
problem_data.h  = @(x, y, iside) zeros ([2, size(x)]); %Dirichlet
%problem_data.g = ... ;%Neumann boundary condition

% Exact solution, to compute the errors
uxex = @(x, y) (2.*exp(x).*(-1 + x).^2.*x.^2.*(-1 + y).*y.*(-1 + 2.*y));
uyex = @(x, y) (-exp(x).*(-1 + x).*x.*(-2 + x.*(3 + x)).*(-1 + y).^2.*y.^2);

problem_data.velex = @(x, y) cat(1, ...
                  reshape (uxex (x,y), [1, size(x)]), ...
                  reshape (uyex (x,y), [1, size(x)]));

problem_data.pressex   = @(x, y) (-424 + 156.*exp(1) + (-1 + y).*y.* ...
              (-456 + exp(x).*(456 + x.^2.*(228 - 5.*(-1 + y).*y) ...
               + 2.*x.*(-228 + (-1 + y).*y) ...
               + 2.*x.^3.*(-36 + (-1 + y).*y) ...
               + x.^4.*(12 + (-1 + y).*y))));

problem_data.gradvelex = @test_stokes_square_graduex;



% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.regularity   = [ 2  2];       % Regularity of the splines
method_data.nsub         = [10 10];       % Number of subdivisions
method_data.nquad        = [ 4  4];       % Points for the Gaussian quadrature rule


% 3) CALL TO THE SOLVER
[geometry, msh, space_v, vel, space_p, press] = ...
                       solve_stokes_2d_bspline_rt (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) COMPARISON WITH EXACT SOLUTION
error_l2_p = sp_l2_error (space_p, msh, press, problem_data.pressex)
[error_h1_v, error_l2_v] = ...
   sp_h1_error (space_v, msh, vel, problem_data.velex, problem_data.gradvelex)


% 4.1) EXPORT TO PARAVIEW
output_file = 'SQUARE_RT_Reg2_Sub10'

vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)'};
sp_to_vtk_2d (press, space_p, geometry, vtk_pts, [output_file '_press'], 'press')
sp_to_vtk_2d (vel,   space_v, geometry, vtk_pts, [output_file '_vel'  ], 'vel')

% 4.3) PLOT IN MATLAB
[eu, F] = sp_eval_2d (vel, space_v, geometry, vtk_pts);
[X,  Y] = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

figure()
subplot(1,2,2)
eu2 = problem_data.velex (X, Y);
quiver (X, Y, squeeze(eu2(1,:,:)), squeeze(eu2(2,:,:)))
axis equal
title('Exact solution')
subplot(1,2,1)
quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
axis equal
title('Computed solution')
