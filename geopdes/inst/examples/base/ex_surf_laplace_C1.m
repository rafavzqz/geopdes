% EX_LAPLACE_LSHAPED_MP: solve the Poisson problem in the multipatch L-shaped domain with a B-spline discretization.
clear all
close all

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  

%Two patch flat planar
% coefs_1(:,1,:)=[0.125 -0.125 0;0.8 0.4 0];
% coefs_1(:,2,:)=[-0.75 0.5 0;0.0 1.0 0];
% 
% coefs_2(:,1,:)=[0.125 -0.125 0;-0.1 -0.9 0];
% coefs_2(:,2,:)=[0.8 0.4 0;0.87 -0.5 0];
% 
% knots{1} = [0 0 1 1]; 
% knots{2} = [0 0 1 1];
% 
% butterf1 = nrbmak(permute(coefs_1,[3,2,1]),knots);
% butterf2 = nrbmak(permute(coefs_2,[3,2,1]),knots);
% 
% nrb(1)=butterf1;
% nrb(2)=butterf2;
% problem_data.geo_name = nrb;


%Three patches flat
% coefs_1(:,3,:)=[0.125 -0.125 1; 0.4625 0.1375 1; 0.8 0.4 1];
% coefs_1(:,2,:)=[-0.3125 0.1875 1; 0.04375 0.44375 1; 0.4 0.7 1];
% coefs_1(:,1,:)=[-0.75 0.5 1;-0.375 0.75 1;  0.0 1.0 1];
% 
% coefs_2(:,3,:)=[0.125 -0.125 1; 0.0125 -0.5125 1; -0.1 -0.9 1];
% coefs_2(:,2,:)=[0.4625 0.1375 1;0.42375 -0.28125 1; 0.385 -0.7 1];
% coefs_2(:,1,:)=[0.8 0.4 1; 0.835 -0.05 1;0.87 -0.5 1];
% 
% coefs_3(:,3,:)=[0.125 -0.125 1; -0.3125 0.1875 1; -0.75 0.5 1];
% coefs_3(:,2,:)=[0.0125 -0.5125 1; -0.39875 -0.25625 1;-0.81 0.0 1];
% coefs_3(:,1,:)=[-0.1 -0.9 1; -0.485 -0.7 1; -0.87 -0.5 1];
% 
% knots{1} = [0 0 0 1 1 1];
% knots{2} = [0 0 0 1 1 1];
% 
% butterf1 = nrbmak(permute(coefs_1,[3,2,1]),knots);
% butterf2 = nrbmak(permute(coefs_2,[3,2,1]),knots);
% butterf3 = nrbmak(permute(coefs_3,[3,2,1]),knots);
% 
% nrb(1)=butterf1;
% nrb(2)=butterf2;
% % nrb(1)=nrbtform (nrb(1), vecrot(pi/3,[1 2 3]));
% % nrb(2)=nrbtform (nrb(2), vecrot(pi/3,[1 2 3]));
% nrb(3)=butterf3;
% problem_data.geo_name = nrb;

% % 3 patch surf
coefs_1(:,3,:)=[0.125 -0.125 0.96875; 0.4625 0.1375 0.95; 0.8 0.4 0.2];
coefs_1(:,2,:)=[-0.3125 0.1875 1.15625; 0.04375 0.44375 1.2625; 0.4 0.7 0.6];
coefs_1(:,1,:)=[-0.75 0.5 0.1875;-0.375 0.75 0.5;  0.0 1.0 0.0];

coefs_2(:,3,:)=[0.125 -0.125 0.96875; 0.0125 -0.5125 0.9; -0.1 -0.9 0.18];
coefs_2(:,2,:)=[0.4625 0.1375 0.95;  0.42375 -0.28125 1.134375; 0.385 -0.7 0.637];
coefs_2(:,1,:)=[0.8 0.4 0.2; 0.835 -0.05 0.504;  0.87 -0.5 -0.0069];

coefs_3(:,3,:)=[0.125 -0.125 0.96875; -0.3125 0.1875 1.15625; -0.75 0.5 0.1875];
coefs_3(:,2,:)=[0.0125 -0.5125 0.9; -0.39875 -0.25625 1.210625; -0.81 0.0 0.5975];
coefs_3(:,1,:)=[-0.1 -0.9 0.18; -0.485 -0.7 0.463; -0.87 -0.5 -0.0069];

knots{1} = [0 0 0 1 1 1];
knots{2} = [0 0 0 1 1 1];

butterf1 = nrbmak(permute(coefs_1,[3,2,1]),knots);
butterf2 = nrbmak(permute(coefs_2,[3,2,1]),knots);
butterf3 = nrbmak(permute(coefs_3,[3,2,1]),knots);

nrb(1)=butterf1;
nrb(2)=butterf2;
nrb(3)=butterf3;
% nrb(1)=nrbtform (nrb(1), vecrot(pi/3,[1 2 3]));
% nrb(2)=nrbtform (nrb(2), vecrot(pi/3,[1 2 3]));
% nrb(3)=nrbtform (nrb(3), vecrot(pi/3,[1 2 3]));
problem_data.geo_name = nrb;


%5 patch surf
% coefs_1(:,1,:)=[0.125 -0.125 0.96875; 0.4375 0.2875 0.99375; 0.75 0.7 -0.0525];
% coefs_1(:,2,:)=[-0.2375 0.3375 1.175; 0.06875 0.59375 1.0075; 0.375 0.85 0.3];
% coefs_1(:,3,:)=[-0.6 0.8 0.0; -0.3 0.9 0.2; 0.0 1.0 0.0];
% 
% coefs_2(:,1,:)=[0.125 -0.125 0.96875; -0.2375 0.3375 1.175; -0.6 0.8 0.0];
% coefs_2(:,2,:)=[-0.3875 -0.1625 1.0875; -0.58125 0.19625 0.88875; -0.775 0.555 0.182];
% coefs_2(:,3,:)=[-0.9 -0.2 0.15; -0.925 0.055 0.207; -0.95 0.31 0.0014];
% 
% coefs_3(:,1,:)=[0.125 -0.125 0.96875; -0.3875 -0.1625 1.0875; -0.9 -0.2 0.15];
% coefs_3(:,2,:)=[0.1125 -0.5458333333333333 0.8666666666666667; -0.29375 -0.5354166666666667 0.9264583333333334; -0.7 -0.525 0.38];
% coefs_3(:,3,:)=[0.1 -0.9666666666666667 0.05555555555555555; -0.2 -0.9083333333333333 0.22833333333333333; -0.5 -0.85 0.0275];
% 
% coefs_4(:,1,:)=[0.125 -0.125 0.96875; 0.1125 -0.5458333333333333 0.8666666666666667; 0.1 -0.9666666666666667 0.05555555555555555];
% coefs_4(:,2,:)=[0.5125 -0.2625 0.8375; 0.40625 -0.5854166666666667 0.6772916666666666; 0.3 -0.9083333333333333 0.12833333333333333];
% coefs_4(:,3,:)=[0.9 -0.4 0.03; 0.7 -0.625 0.21; 0.5 -0.85 0.0275];
% 
% coefs_5(:,1,:)=[0.125 -0.125 0.96875; 0.5125 -0.2625 0.8375; 0.9 -0.4 0.03];
% coefs_5(:,2,:)=[0.4375 0.2875 0.99375; 0.68125 0.12125 0.7625; 0.925 -0.045 0.269];
% coefs_5(:,3,:)=[0.75 0.7 -0.0525; 0.85 0.505 0.0705; 0.95 0.31 0.0014];
% 
% knots{1} = [0 0 0 1 1 1];
% knots{2} = [0 0 0 1 1 1];
% 
% butterf1 = nrbmak(permute(coefs_1,[3,2,1]),knots);
% butterf2 = nrbmak(permute(coefs_2,[3,2,1]),knots);
% butterf3 = nrbmak(permute(coefs_3,[3,2,1]),knots);
% butterf4 = nrbmak(permute(coefs_4,[3,2,1]),knots);
% butterf5 = nrbmak(permute(coefs_5,[3,2,1]),knots);
% 
% nrb(1)=butterf1;
% nrb(2)=butterf2;
% nrb(3)=butterf3;
% nrb(4)=butterf4;
% nrb(5)=butterf5;
% problem_data.geo_name = nrb;


% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [];
problem_data.weak_drchlt_sides = [1 2 3 4 5 6];

% Physical parameters
problem_data.c_diff  = @(x, y, z) ones(size(x));

% % Source and boundary terms
% problem_data.f = @(x, y, z) exp(x).*((x.^2 + y.^2 - 1).*sin(x.*y) - 2*y.*cos(x.*y));
% problem_data.g = @(x, y, z, ind) test_Lshaped_mp_g_nmnn (x, y, ind);
% problem_data.h = @(x, y, z, ind) exp(x) .* sin (x.*y);

% % Exact solution (optional)
% problem_data.uex     = @(x, y, z) exp(x) .* sin (x.*y);
% problem_data.graduex = @(x, y, z) cat (1, ...
%                reshape (exp(x).*(sin(x.*y) + y.*cos(x.*y)), [1, size(x)]), ...
%                reshape (exp(x).*x.*cos(x.*y), [1, size(x)]));

 
% Surf
% problem_data.f = @(x, y, z) zeros(size(x));
% problem_data.h = @(x, y, z, ind) exp(y) .* sin (x);
% problem_data.uex     = @(x, y, z) exp (y) .* sin (x);
% problem_data.graduex = @(x, y, z) cat (1, ...
%                        reshape (exp(y).*cos(x), [1, size(x)]), ...
%                        reshape (exp(y).*sin(x), [1, size(x)]), ...
%                         zeros([1, size(x)]));
                    
                    
% problem_data.f = @(x, y, z) -(1 ./ ( (1 + 4 * x.^2 + 4 * y.^2).^2) ) .* (-4 * cos(9 * x) .* ( (65 + 392 * x.^4 + 422 * y.^2 + 648 * y.^4 + ...
%                             2 * x.^2 .* (179 + 520 * y.^2) ) .* cos(1 - 7 * y) + 28 * y .* (1 + 2 * x.^2 + 2 * y.^2) .* sin(1 - 7 * y) ) + ...
%                             144 * x .* sin(9 * x) .* ( (1 + 2 * x.^2 + 2 * y.^2) .* cos(1 - 7 * y) + 7 * y .* (1 + 4 * x.^2 + 4 * y.^2) .* sin(1 - 7 * y) ) );
% problem_data.h = @(x, y, z, ind) cos(9*x) .* cos(1 - 7*y);
% problem_data.uex     = @(x, y, z) cos(9*x) .* cos(1 - 7*y);
% problem_data.graduex = @(x, y, z) cat (1, ...
%                        reshape (-( (2 * (9 * (1 + 4 * y.^2) .* cos(1 - 7 * y) .* sin(9 * x) + 28 * x .* y .* cos(9 * x) .* sin(1 - 7 * y) ) ) ./ (1 + 4 * x.^2 + 4 * y.^2)), [1, size(x)]), ...
%                        reshape ((72 * x .* y .* cos(1 - 7 * y) .* sin(9 * x) + 14 * (1 + 4 * x.^2) .* cos(9 * x) .* sin(1 - 7 * y) ) ./ (1 + 4 * x.^2 + 4 * y.^2), [1, size(x)]), ...
%                        reshape((4 * (9 * x .* cos(1 - 7 * y) .* sin(9 * x) - 7 * y .* cos(9 * x) .* sin(1 - 7 * y) ) ) ./ (1 + 4 * x.^2 + 4 * y.^2), [1, size(x)]));                 
                    


% problem_data.f = @(x, y, z) ( (4 * (1 + 2 * x.^2 + 2 * y.^2) ) ./ (1 + 4 * x.^2 + 4 * y.^2).^2);
% problem_data.h = @(x, y, z, ind) z;
% problem_data.uex     = @(x, y, z) z;
% problem_data.graduex = @(x, y, z) cat (1, ...
%                        reshape (-( (2 * x) ./ (1 + 4 * x.^2 + 4 * y.^2) ), [1, size(x)]), ...
%                        reshape (-( (2 * y) ./ (1 + 4 * x.^2 + 4 * y.^2) ), [1, size(x)]), ...
%                        reshape ( (4 * (x.^2 + y.^2) ) ./ (1 + 4 * x.^2 + 4 * y.^2), [1, size(x)]));                             
                   
                  
problem_data.f = @(x, y, z) (2./((1 + 4*x.^2 + 4*y.^2).^2)).*(cos(3*x).*((5 + 8*x.^4 + 38*y.^2 + 72*y.^4 + x.^2 .*(22 + 80*y.^2)).*cos(y) - 4*y.*(1 + 2*x.^2 + 2*y.^2).*sin(y)) + 12*x.*sin(3*x).*(-(1 + 2*x.^2 + 2*y.^2).*cos(y) + y.*(1 + 4*x.^2 + 4*y.^2).*sin(y)));
problem_data.h = @(x, y, z, ind) cos(3*x) .* cos(y);
problem_data.uex     = @(x, y, z) cos(3*x) .* cos(y);
problem_data.graduex = @(x, y, z) cat (1, ...
                       reshape ((-3*(1 + 4*y.^2).*cos(y).*sin(3*x) + 4*x.*y.*cos(3*x).*sin(y))./(1 + 4*x.^2 + 4*y.^2), [1, size(x)]), ...
                       reshape ((12*x.*y.*cos(y).*sin(3*x) - (1 + 4*x.^2).*cos(3*x).*sin(y))./(1 + 4*x.^2 + 4*y.^2), [1, size(x)]), ...
                       reshape ((2*(3*x.*cos(y).*sin(3*x) + y.*cos(3*x).*sin(y)))./(1 + 4*x.^2 + 4*y.^2), [1, size(x)]));                
                   
                   
% problem_data.f = @(x, y) zeros (size(x));
% problem_data.h = @(x, y, ind) exp(y) .* sin (x);
% problem_data.uex     = @(x, y) exp (y) .* sin (x);
% problem_data.graduex = @(x, y) cat (1, ...
%                        reshape (exp(y).*cos(x), [1, size(x)]), ...
%                        reshape (exp(y).*sin(x), [1, size(x)]));
                    
             
% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [5 5];       % Degree of the splines
method_data.regularity = [1 1];       % Regularity of the splines

method_data.nquad      = method_data.degree+1;       % Points for the Gaussian quadrature rule
method_data.Cpen = 10 * (min(method_data.degree) + 1);

h=[];

for i= 1:4
    fprintf ('Loop %d \n', i);
    method_data.nsub       = [2^(i+1) 2^(i+1)] ;      % Number of subdivisions
    % 3) CALL TO THE SOLVER
    [geometry, msh, space, u] = mp_solve_laplace_C1 (problem_data, method_data);

    h = [h 1/sqrt(msh.nel_per_patch(1))]
    % % 4) POST-PROCESSING
    % % EXPORT TO PARAVIEW
    % output_file = 'Lshaped_mp_BSP_Deg3_Reg2_Sub9';
    % 
    vtk_pts = {linspace(0, 1, 50), linspace(0, 1, 50)};
    % fprintf ('The result is saved in the file %s.pvd \n \n', output_file);
    % sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')
    % 

    % COMPARISON WITH THE EXACT SOLUTION
    [error_h1(i), error_l2(i)] = sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex)
end

figure
sp_plot_solution (u, space, geometry, vtk_pts)

figure
loglog(h, error_l2,'-o')
hold on
loglog(h, error_h1,'-*')
loglog(h, 0.0001*h.^6,'-x')
loglog(h, 0.001*h.^5,'-x')
legend("L^2 error", "H^1 error", "h^6", "h^5",'Location','southeast');

