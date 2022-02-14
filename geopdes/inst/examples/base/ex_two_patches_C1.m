% EX_LAPLACE_LSHAPED_MP: solve the Poisson problem in the multipatch L-shaped domain with a B-spline discretization.
clear all
close all
% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
% nrb(1) = nrb4surf ([0 0], [-1 0], [0 1], [-1 1]);
% nrb(2) = nrb4surf ([0 0], [1 0], [0 1], [1 1]);

% p1 = [-4 -1/2]; p2 = [0 0]; p3 = [-10/3 16/5]; p4 = [0 3];
% nrb(1) = nrb4surf (p2, p1, p4, p3);
% 
% p1 = [0 0]; p2 = [8/3 -2/5]; p3 = [0 3]; p4 = [10/3 23/7];
% nrb(2) = nrb4surf (p1, p2, p3, p4); problem_data.geo_id = 1;
% problem_data.geo_name = nrb;

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

% 3 patch surf
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
problem_data.geo_name = nrb;


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

 
%Surf
problem_data.f = @(x, y, z) zeros (size(x));
problem_data.h = @(x, y, z, ind) exp(y) .* sin (x);
problem_data.uex     = @(x, y, z) exp (y) .* sin (x);
problem_data.graduex = @(x, y, z) cat (1, ...
                       reshape (exp(y).*cos(x), [1, size(x)]), ...
                       reshape (exp(y).*sin(x), [1, size(x)]), ...
                        zeros([1, size(x)]));
                    
% problem_data.f = @(x, y) zeros (size(x));
% problem_data.h = @(x, y, ind) exp(y) .* sin (x);
% problem_data.uex     = @(x, y) exp (y) .* sin (x);
% problem_data.graduex = @(x, y) cat (1, ...
%                        reshape (exp(y).*cos(x), [1, size(x)]), ...
%                        reshape (exp(y).*sin(x), [1, size(x)]));
                    
             
% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3];       % Degree of the splines
method_data.regularity = [1 1];       % Regularity of the splines
method_data.nsub       = [16 16];       % Number of subdivisions
method_data.nquad      = method_data.degree+1;       % Points for the Gaussian quadrature rule
method_data.Cpen = 10 * (min(method_data.degree) + 1);

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = mp_solve_laplace_C1 (problem_data, method_data);

% % 4) POST-PROCESSING
% % EXPORT TO PARAVIEW
% output_file = 'Lshaped_mp_BSP_Deg3_Reg2_Sub9';
% 
vtk_pts = {linspace(0, 1, 50), linspace(0, 1, 50)};
% fprintf ('The result is saved in the file %s.pvd \n \n', output_file);
% sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')
% 
figure
sp_plot_solution (u, space, geometry, vtk_pts)

% COMPARISON WITH THE EXACT SOLUTION
[error_h1, error_l2] = sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex)

