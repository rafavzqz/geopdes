% EX_STOKES_ANNULUS_RT: solve the Stokes problem in one quarter of an annulus with generalized Raviart-Thomas elements.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'annulus.mat';

% Type of boundary conditions for each side of the domain
problem_data.drchlt_sides = 1:4;
problem_data.nmnn_sides = [];

% Physical parameters
problem_data.viscosity = @(x, y) ones (size (x));

% Force term
fx = @(x, y) (4*(20*(98-57*x.^2).*y.^7+42*x.*(34*x.^2-75).*y.^6 ...
                 -3*(430*x.^4-1480*x.^2+1023).*y.^5 ...
                 +15*x.*(70*x.^4-310*x.^2+297).*y.^4 ...
                 +3*x.*(x.^4-5*x.^2+4).^2 ...
                 +2*(-290*x.^6+1500*x.^4-2079*x.^2+680).*y.^3 ...
                 +6*x.*(38*x.^6-255*x.^4+495*x.^2-260).*y.^2 ...
                 +(-75*x.^8+520*x.^6-1089*x.^4+720*x.^2-112).*y ...
                 +603*x.*y.^8-355*y.^9)-y./(x.^2+y.^2));

fy = @(x, y) (x./(x.^2+y.^2)+4*(5*x.^9-27*x.^8.*y+20*x.^7.*(15*y.^2-2) ...
                                +14*x.^6.*y.*(15-38*y.^2)+3*x.^5.*(290*y.^4-520*y.^2+33) ...
                                -15*x.^4.*y.*(70*y.^4-170*y.^2+33) ...
                                +x.^3.*(860*y.^6-3000*y.^4+2178*y.^2-80) ...
                                +18*x.^2.*y.*(-34*y.^6+155*y.^4-165*y.^2+20) ...
                                +x.*(285*y.^8-1480*y.^6+2079*y.^4-720*y.^2+16) ...
                                -67*y.^9+450*y.^7-891*y.^5+520*y.^3-48*y));
problem_data.f  = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));

% Boundary terms
problem_data.h  = @(x, y, iside) zeros ([2, size(x)]); %Dirichlet

% Exact solution, to compute the errors
uxex = @(x, y) (2*(x-y).*y.*(-4+x.^2+y.^2).*(-1+x.^2+y.^2).* ...
                (x.^5-8*y-2*x.^4.*y+20*y.^3-6*y.^5+2*x.^2.*y.*(5-4*y.^2)+x.^3.*(-5 + 6*y.^2)+ ...
                 x.*(4+5*y.^2.*(-3+y.^2))));

uyex = @(x, y) (-2*(x - y).*y.^2.*(-4 + x.^2 + y.^2).*(-1 + x.^2 + y.^2).* ...
                (4 + 5*x.^2.*(-3 + x.^2) + 2*x.*(5 - 2*x.^2).*y + ...
                 (-5 + 6*x.^2).*y.^2 - 4*x.*y.^3 + y.^4));
problem_data.velex = @(x, y) cat(1, ...
                  reshape (uxex (x,y), [1, size(x)]), ...
                  reshape (uyex (x,y), [1, size(x)]));

problem_data.gradvelex = @test_stokes_annulus_graduex;

problem_data.pressex = @(x, y) (-(pi/8) + atan (y./x));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.element_name = 'rt';     % Element type for discretization
method_data.degree       = [ 3  3];  % Degree of the splines
method_data.regularity   = [ 2  2];  % Regularity of the splines
method_data.nsub         = [10 10];  % Number of subdivisions
method_data.nquad        = [ 5  5];  % Points for the Gaussian quadrature rule

% Penalization parameter for Nitsche's method
factor = 10;
method_data.Cpen = factor*(min(method_data.degree)+1);

% 3) CALL TO THE SOLVER
[geometry, msh, space_v, vel, space_p, press] = ...
                       solve_stokes_2d (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) COMPARISON WITH EXACT SOLUTION
error_l2_p = sp_l2_error (space_p, msh, press, problem_data.pressex)
[error_h1_v, error_l2_v] = ...
   sp_h1_error (space_v, msh, vel, problem_data.velex, problem_data.gradvelex)


% 4.2) EXPORT TO PARAVIEW
output_file = 'ANNULUS_RT_Deg3_Reg2_Sub10';

fprintf ('The result is saved in the files %s \n and %s \n \n', ...
           [output_file '_vel'], [output_file '_press']);
vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
sp_to_vtk (press, space_p, geometry, vtk_pts, [output_file '_press'], 'press')
sp_to_vtk (vel,   space_v, geometry, vtk_pts, [output_file '_vel'  ], 'vel')

% 4.3) PLOT IN MATLAB
[eu, F]  = sp_eval (vel, space_v, geometry, vtk_pts);
[X,  Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

figure()
subplot (1,2,2)
eu2 = problem_data.velex (X, Y);
quiver (X, Y, squeeze(eu2(1,:,:)), squeeze(eu2(2,:,:)))
axis equal
title('Exact solution')
subplot (1,2,1)
quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
axis equal
title('Computed solution')

[div, F] = sp_eval (vel, space_v, geometry, vtk_pts, 'divergence');
figure()
surf (X, Y, div)
view(2)
axis equal
title('Computed divergence')

%!demo
%! ex_stokes_annulus_rt