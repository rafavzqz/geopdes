% EX_NAVIERSTOKES_SQUARE: test convergence of Navier-Stokes solver with an
% analytical solution defined on a square geometry.
clear problem_data
nrb = nrbsquare([0 0],1,1);
problem_data.geo_name = nrb;

% boundary conditions
problem_data.drchlt_sides = [1 2];
problem_data.nmnn_sides = [3 4];

% Physical parameters
mu = 1;
problem_data.viscosity = @(x, y) mu * ones (size (x));

% Force term
u1ex = @(x,y) cos(y*pi);
u2ex = @(x,y) x.*(x-1);
pex = @(x,y)  sin(pi*y);

u1exdx = @(x,y) 0*x;
u1exdy = @(x,y) -pi*sin(y*pi);
u1exdxdx = @(x,y) 0*x;
u1exdydy = @(x,y) -pi^2*cos(y*pi);

u2exdx = @(x,y) 2*x -x.^0;
u2exdy = @(x,y) 0*x;
u2exdxdx = @(x,y) 0*x;
u2exdydy = @(x,y) 2*x.^0;

pexdx = @(x,y) 0*x;
pexdy = @(x,y) pi*cos(pi*y);

fx = @(x, y) -mu*(u1exdxdx(x,y)+u1exdydy(x,y)) + (u1ex(x,y).*u1exdx(x,y) + u2ex(x,y).*u1exdy(x,y)) + pexdx(x,y);
fy = @(x, y) -mu*(u2exdxdx(x,y)+u2exdydy(x,y)) + (u1ex(x,y).*u2exdx(x,y) + u2ex(x,y).*u2exdy(x,y)) + pexdy(x,y);

problem_data.f  = @(x, y) cat(1, ...
    reshape (fx (x,y), [1, size(x)]), ...
    reshape (fy (x,y), [1, size(x)]));

% Boundary terms
problem_data.h  = @(x, y, iside) cat(1,reshape(u1ex(x,y),[1,size(x)]), ...
                                       reshape(u2ex(x,y),[1,size(x)]));
problem_data.g  = @(x, y, iside) cat(1,zeros([1,size(x)]),zeros([1,size(x)]));

% Exact solution, to compute the errors
problem_data.velex = @(x, y) cat(1, ...
    reshape (u1ex (x,y), [1, size(x)]), ...
    reshape (u2ex (x,y), [1, size(x)]));

problem_data.gradvelex = @test_navierstokes_square_graduex;

problem_data.pressex = @(x, y) sin(pi*y);

deg = 2;

for ii = 1:4
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Refinement iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n', ii);
    clear method_data
    method_data.element_name = 'RT';     % Element type for discretization
    method_data.degree       = [deg deg];  % Degree of the splines (pressure space)
    method_data.regularity   = [deg-1 deg-1];  % Regularity of the splines (pressure space)
    method_data.nsub         = [2^ii  2^ii];  % Number of subdivisions
    method_data.nquad        = [deg+1 deg+1];  % Points for the Gaussian quadrature rule
    
    % Penalization parameter for Nitsche's method
    factor = 10;
    method_data.Cpen = factor*(min(method_data.degree)+1);
    
    solver_data.name = 'newton';
    solver_data.tol = 1e-8;
    solver_data.nmax = 20;
    
    [geometry, msh, space_v, vel, space_p, press] = ...
        solve_navier_stokes (problem_data, method_data,solver_data);
    
    error_l2_p(ii) = sp_l2_error (space_p, msh, press, problem_data.pressex);
    [error_h1_v(ii), error_l2_v(ii)] = ...
        sp_h1_error (space_v, msh, vel, problem_data.velex, problem_data.gradvelex);
    
    msh_prc = msh_precompute(msh);
    h(ii) = max(msh_prc.element_size);
end

figure
loglog(h, error_h1_v, 'b-o', 'Markersize', 8, 'Linewidth', 2);
hold on
grid on
loglog(h, error_l2_p, 'r-x', 'Markersize', 8, 'Linewidth', 2);
loglog(h, 0.5*h.^(deg+1), '-k', 'LineWidth', 2);
xlabel('Mesh size');
ylabel('Error');
slope = strcat('h^', num2str(deg+1));
if (strcmpi(method_data.element_name,'RT'))
  loglog(h, h.^deg, '--k', 'LineWidth', 2);
  slope_v = strcat('h^', num2str(deg));
  legend('Velocity, H1 error', 'Pressure, L2 error', slope, slope_v);
else
  legend('Velocity, H1 error', 'Pressure, L2 error', slope);
end
