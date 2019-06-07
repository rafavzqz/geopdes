% TEST_CONVERGENCE: test convergence of Navier-Stokes solver with an
% analytical solution defined on a square geometry.
clc
clear all
close all

clear problem_data
nurb = nrbsquare([0 0],1,1);
problem_data.geo_name = nurb;

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

problem_data.gradvelex = @test_navier_stokes_convergence_graduex;

problem_data.pressex = @(x, y) sin(pi*y);

deg = 2;

for i = 1:4
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Refinement iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n', i);
    clear method_data
    method_data.element_name = 'rt';     % Element type for discretization
    method_data.degree       = [deg deg];  % Degree of the splines (pressure space)
    method_data.regularity   = [deg-1 deg-1];  % Regularity of the splines (pressure space)
    method_data.nsub         = [3*2^i  3*2^i];  % Number of subdivisions
    method_data.nquad        = [deg+1 deg+1];  % Points for the Gaussian quadrature rule
    
    % Penalization parameter for Nitsche's method
    factor = 10;
    method_data.Cpen = factor*(min(method_data.degree)+1);
    
    solver_data.name = 'newton';
    solver_data.tol = 1e-8;
    solver_data.nmax = 20;
    
    [geometry, msh, space_v, vel, space_p, press] = ...
        solve_navier_stokes (problem_data, method_data,solver_data);
    
    error_l2_p(i) = sp_l2_error (space_p, msh, press, problem_data.pressex);
    [error_h1_v(i), error_l2_v(i)] = ...
        sp_h1_error (space_v, msh, vel, problem_data.velex, problem_data.gradvelex);
    
    msh_prc = msh_precompute(msh);
    h(i) = max(msh_prc.element_size);
end


if(strcmp(method_data.element_name,'rt'))
    figure
    loglog(h,error_h1_v,'.-','Markersize',10,'Linewidth',1);
    title('H1-error on velocity Raviart-Thomas element');
    hold on
    loglog(h,h.^deg,'--k');
    grid on
    xlabel('h');
    ylabel('error');
    slope = strcat('h^', num2str(deg));
    legend('error',slope);

    figure
    loglog(h,error_l2_p,'.-','Markersize',10,'Linewidth',1);
    title('L2-error on pressure Raviart-Thomas element');
    hold on
    loglog(h,h.^(deg+1)*1e-2,'--k');
    grid on
    xlabel('h');
    ylabel('error');
    slope = strcat('h^', num2str(deg+1));
    legend('error',slope);
elseif(strcmp(method_data.element_name,'th'))
    figure
    loglog(h,error_h1_v,'.-','Markersize',10,'Linewidth',1);
    title('H1-error on velocity Taylor-Hood element');
    hold on
    loglog(h,h.^(deg+1),'--k');
    grid on
    xlabel('h');
    ylabel('error');
    slope = strcat('h^', num2str(deg));
    legend('error',slope);

    figure
    loglog(h,error_l2_p,'.-','Markersize',10,'Linewidth',1);
    title('L2-error on pressure Taylor-Hood element');
    hold on
    loglog(h,h.^(deg+1)*1e-2,'--k');
    grid on
    xlabel('h');
    ylabel('error');
    slope = strcat('h^', num2str(deg+1));
    legend('error',slope);

elseif(strcmp(method_data.element_name,'NDL'))
    figure
    loglog(h,error_h1_v,'.-','Markersize',10,'Linewidth',1);
    title('H1-error on velocity Nedelec element');
    hold on
    loglog(h,h.^(deg+1),'--k');
    grid on
    xlabel('h');
    ylabel('error');
    slope = strcat('h^', num2str(deg));
    legend('error',slope);

    figure
    loglog(h,error_l2_p,'.-','Markersize',10,'Linewidth',1);
    title('L2-error on pressure Nedelec element');
    hold on
    loglog(h,h.^(deg+1)*1e-2,'--k');
    grid on
    xlabel('h');
    ylabel('error');
    slope = strcat('h^', num2str(deg+1));
    legend('error',slope);
elseif(strcmp(method_data.element_name,'SG'))
    figure
    loglog(h,error_h1_v,'.-','Markersize',10,'Linewidth',1);
    title('H1-error on velocity Sub-grid element');
    hold on
    loglog(h,h.^(deg+1)*6e-1,'--k');
    grid on
    xlabel('h');
    ylabel('error');
    slope = strcat('h^', num2str(deg));
    legend('error',slope);

    figure
    loglog(h,error_l2_p,'.-','Markersize',10,'Linewidth',1);
    title('L2-error on pressure Sub-grid element');
    hold on
    loglog(h,h.^(deg+1)*3e-1,'--k');
    grid on
    xlabel('h');
    ylabel('error');
    slope = strcat('h^', num2str(deg+1));
    legend('error',slope);
   
else
    error('Unknown element type')
end


function gu = test_navier_stokes_convergence_graduex (x, y)
  uxx = @(x, y) 0*x;
  uxy = @(x, y) -pi*sin(y*pi);
  uyx = @(x, y) 2*x - x.^0;
  uyy = @(x, y) 0*x;
  
  
  gu = zeros (2, 2, size(x,1), size(x,2));
  gu(1, 1, :, :) = reshape (uxx (x,y), [1, 1, size(x)]);
  gu(2, 1, :, :) = reshape (uyx (x,y), [1, 1, size(x)]);
  gu(1, 2, :, :) = reshape (uxy (x,y), [1, 1, size(x)]);
  gu(2, 2, :, :) = reshape (uyy (x,y), [1, 1, size(x)]);
end
